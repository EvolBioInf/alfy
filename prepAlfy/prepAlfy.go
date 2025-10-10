package main

import (
	"bytes"
	"flag"
	"fmt"
	"github.com/evolbioinf/alfy/mat"
	"github.com/evolbioinf/alfy/util"
	"github.com/evolbioinf/clio"
	"github.com/evolbioinf/esa"
	"github.com/evolbioinf/fasta"
	"log"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"slices"
	"sort"
	"strconv"
	"strings"
	"sync"
)

type Matches struct {
	matchLengths [][]int
	subjectID    [][][]int
	id           int
}

func readDir(dir string) map[string]bool {
	dirEntries, err := os.ReadDir(dir)
	util.Check(err)
	names := make(map[string]bool)
	for _, dirEntry := range dirEntries {
		if dirEntry.IsDir() {
			p := dir + "/" + dirEntry.Name()
			fmt.Fprintf(os.Stderr,
				"Skipping subdirectory %s\n", p)
			continue
		}
		ext := filepath.Ext(dirEntry.Name())
		if ext != ".fasta" && ext != ".fna" &&
			ext != ".ffn" && ext != ".fnr" &&
			ext != ".fa" {
			m := "%s doesnt have the extension of a " +
				"nucleotide FASTA file." +
				"Skipping it \n"
			p := dir + "/" + dirEntry.Name()
			fmt.Fprintf(os.Stderr, m, p)
			continue
		}
		names[dirEntry.Name()] = true
	}
	return names
}
func main() {
	util.PrepareErrorMessages("prepAlfy")
	optV := flag.Bool("v", false, "version")
	optQ := flag.String("q", "", "query directory")
	optS := flag.String("s", "", "subject directory")
	ncpu := runtime.NumCPU()
	optT := flag.Int("t", ncpu, "number of threads")
	optN := flag.Bool("n", false, "print subject names "+
		"(default identifiers)")
	optR := flag.Bool("r", false, "raw output "+
		"(default intervals)")
	optSS := flag.Bool("S", false, "step (default skip)")
	u := "prepAlfy -q <queryDir> -s <subjectDir>"
	p := "Prepare alfy input"
	e := "prepAlfy -q queries/ -s subjects/"
	clio.Usage(u, p, e)
	flag.Parse()
	if *optV {
		util.Version("prepAlfy")
	}
	if *optQ == "" {
		m := "please provide a directory " +
			"of query sequences"
		fmt.Fprintf(os.Stderr, "%s\n", m)
		os.Exit(1)
	}
	if *optS == "" {
		m := "please provide a directory " +
			"of subject sequences"
		fmt.Fprintf(os.Stderr, "%s\n", m)
		os.Exit(1)
	}
	if *optT <= 0 {
		log.Fatalf("Can't set %d threads.", *optT)
	}
	queries := readDir(*optQ)
	if len(queries) == 0 {
		fmt.Fprintf(os.Stderr, "%s is empty\n", *optQ)
		os.Exit(1)
	}
	subjects := readDir(*optS)
	if len(subjects) == 0 {
		fmt.Fprintf(os.Stderr, "%s is empty\n", *optS)
		os.Exit(1)
	}
	var queryNames, subjectNames []string
	var subjectIDs []int
	for query := range queries {
		queryNames = append(queryNames, query)
	}
	sort.Strings(queryNames)
	for subject := range subjects {
		subjectNames = append(subjectNames, subject)
	}
	sort.Strings(subjectNames)
	for i, _ := range subjectNames {
		subjectIDs = append(subjectIDs, i)
	}
	for _, query := range queryNames {
		if subjects[query] {
			m := "Found %s%s and %s%s." +
				"Please ensure queries and " +
				"subjects do not " +
				"overlap."
			fmt.Fprintf(os.Stderr, m, *optQ, query,
				*optS, query)
			os.Exit(1)
		}
	}
	msl := -1
	for _, subject := range subjectNames {
		p := *optS + "/" + subject
		f, err := os.Open(p)
		util.Check(err)
		sc := fasta.NewScanner(f)
		for sc.ScanSequence() {
			l := len(sc.Sequence().Data())
			if l > msl {
				msl = l
			}
		}
		f.Close()
	}
	sID := make(map[int]string)
	esas := make([]*esa.Esa, len(subjectNames))
	for i, subject := range subjectNames {
		sID[i] = subject
		p := *optS + "/" + subject
		f, err := os.Open(p)
		util.Check(err)
		sc := fasta.NewScanner(f)
		for sc.ScanSequence() {
			s := sc.Sequence()
			d := s.Data()
			h := []byte(s.Header())
			for len(d) < msl && sc.ScanSequence() {
				s = sc.Sequence()
				h = append(h, '|')
				h = append(h, []byte(s.Header())...)
				d = append(d, s.Data()...)
			}
			d = bytes.ToUpper(d)
			esas[i] = esa.MakeEsa(d)
		}
	}
	for _, query := range queryNames {
		p := *optQ + "/" + query
		f, err := os.Open(p)
		util.Check(err)
		var querySeqs, revQuerySeqs []*fasta.Sequence
		sc := fasta.NewScanner(f)
		for sc.ScanSequence() {
			s := sc.Sequence()
			querySeqs = append(querySeqs, s)
			s = fasta.NewSequence(s.Header(), s.Data())
			s.ReverseComplement()
			revQuerySeqs = append(revQuerySeqs, s)
		}
		f.Close()
		for i, querySeq := range querySeqs {
			h := querySeq.Header()
			d := bytes.ToUpper(querySeq.Data())
			querySeqs[i] = fasta.NewSequence(h, d)
			h = revQuerySeqs[i].Header()
			d = bytes.ToUpper(revQuerySeqs[i].Data())
			revQuerySeqs[i] = fasta.NewSequence(h, d)
		}
		EsaSets := make([][]*esa.Esa, 0)
		subjectIDsets := make([][]int, 0)
		n := len(subjectNames)
		length := int(math.Ceil(float64(n) / float64(*optT)))
		start := 0
		end := length
		for start < n {
			EsaSets = append(EsaSets,
				esas[start:end])
			subjectIDsets = append(subjectIDsets,
				subjectIDs[start:end])
			start = end
			end += length
			if end > n {
				end = n
			}
		}
		matchesSets := make(chan Matches)
		var wg sync.WaitGroup
		for i, Esas := range EsaSets {
			subjectIDs := subjectIDsets[i]
			wg.Add(1)
			go func(Esas []*esa.Esa, subjectIDs []int) {
				defer wg.Done()
				var matches Matches
				matches.id = i
				n := len(querySeqs)
				matches.matchLengths = make([][]int, n)
				matches.subjectID = make([][][]int, n)
				for i, querySeq := range querySeqs {
					m := len(querySeq.Data())
					matches.matchLengths[i] = make([]int, m)
					matches.subjectID[i] = make([][]int, m)
					for j := range matches.subjectID[i] {
						matches.subjectID[i][j] = make([]int, 1)
					}
				}
				for i, querySeq := range querySeqs {
					f := querySeq.Data()
					rev := false
					for j, esa := range Esas {
						mat.GetMatchLengths(f, esa, subjectIDs[j],
							matches.matchLengths[i],
							matches.subjectID[i],
							rev,
							(*optSS))
					}
					r := revQuerySeqs[i].Data()
					rev = true
					for j, esa := range Esas {
						mat.GetMatchLengths(r, esa, subjectIDs[j],
							matches.matchLengths[i],
							matches.subjectID[i],
							rev,
							(*optSS))
					}
				}
				matchesSets <- matches
			}(Esas, subjectIDs)
		}
		go func() {
			wg.Wait()
			close(matchesSets)
		}()
		matchLengths := make([][]int, 0)
		subjectIDs := make([][][]int, 0)
		for i, qs := range querySeqs {
			n := len(qs.Data())
			ml := make([]int, n)
			matchLengths = append(matchLengths, ml)
			mls := make([][]int, n)
			subjectIDs = append(subjectIDs, mls)
			for j := range subjectIDs[i] {
				subjectIDs[i][j] = append(subjectIDs[i][j], 0)
			}
		}
		msSlice := []Matches{}
		for matchesSet := range matchesSets {
			msSlice = append(msSlice, matchesSet)
		}
		slices.SortFunc(msSlice, func(a, b Matches) int {
			if a.id < b.id {
				return -1
			} else if a.id == b.id {
				return 0
			}
			return 1
		})
		for _, match := range msSlice {
			for i := range match.matchLengths {
				mat.Update(matchLengths[i],
					match.matchLengths[i],
					subjectIDs[i], match.subjectID[i])
			}
		}
		for i := range matchLengths {
			mat.Interpolate(matchLengths[i],
				subjectIDs[i])
		}
		e := strings.LastIndex(query, ".")
		query = query[:e]
		fmt.Printf("#%s\n", query)
		for i, querySeq := range querySeqs {
			ml := matchLengths[i]
			ids := subjectIDs[i]
			fmt.Printf(">%s", querySeq.Header())
			if !*optN {
				marked := make(map[int]bool)
				seen := []int{}
				for _, id := range ids {
					for _, p := range id {
						marked[p] = true
					}
				}
				for i, id := range marked {
					if id {
						seen = append(seen, i)
					}
				}
				slices.Sort(seen)
				for _, id := range seen {
					fmt.Printf(" %d=%s", id+1, subjectNames[id])
				}
			}
			fmt.Printf("\n")
			if *optR {
				for j := 0; j < len(ml); j++ {
					fmt.Println(j, ml[j], ids[j])
				}
			} else {
				start := make([]int, len(ml))
				end := make([]int, len(ml))
				mat := make([]int, len(ml))
				start[0] = 0
				mat[0] = ml[0]
				row := 0
				for j := 1; j < len(ml); j++ {
					prev := ml[j-1]
					curr := ml[j]
					pId := ids[j-1]
					slices.Sort(pId)
					slices.Sort(ids[j])
					str := make([]string, len(pId))
					if *optN {
						for i, val := range pId {
							name := strconv.Itoa(val + 1)
							name = subjectNames[val]
							e = strings.LastIndex(name, ".")
							name = name[:e]
							str[i] = fmt.Sprintf("%v", name)
						}
					} else {
						for i, val := range pId {
							str[i] = fmt.Sprintf("%d", val+1)
						}
					}
					if curr > prev ||
						!slices.Equal(pId, ids[j]) {
						end[row] = j - 1
						fmt.Printf("%d\t%d\t%d\t%v\n",
							start[row], end[row],
							mat[row], strings.Join(str, ","))
						row++
						start[row] = j
						mat[row] = ml[j]
					}
					if j == len(ml)-1 {
						fmt.Printf("%d\t%d\t%d\t%v\n",
							start[row], len(ml)-1,
							mat[row], strings.Join(str, ","))
					}
				}
			}
		}
	}
}
