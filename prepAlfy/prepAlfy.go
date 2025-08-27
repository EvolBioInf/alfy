package main

import (
	"bufio"
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
	"sort"
	"strings"
	"sync"
)

type Matches struct {
	matchLengths [][]int
	subjectID    [][]int
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
	optT := flag.Int("T", ncpu, "number of threads")
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
	sID := make(map[int]string)
	for i, subject := range subjectNames {
		sID[i] = subject
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
	for query, _ := range queries {
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
		subjectNameSets := make([][]string, 0)
		subjectIDsets := make([][]int, 0)
		n := len(subjectNames)
		length := int(math.Ceil(float64(n) / float64(*optT)))
		start := 0
		end := length
		for start < n {
			subjectNameSets = append(subjectNameSets,
				subjectNames[start:end])
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
		for i, subjectNames := range subjectNameSets {
			subjectIDs := subjectIDsets[i]
			wg.Add(1)
			var matches Matches
			go func(subjectNames []string, subjectIDs []int) {
				defer wg.Done()
				for _, querySeq := range querySeqs {
					n := len(querySeq.Data())
					lengths := make([]int, n)
					matches.matchLengths =
						append(matches.matchLengths, lengths)
					n = len(querySeq.Data())
					lengths = make([]int, n)
					matches.subjectID =
						append(matches.subjectID, lengths)
					for j, subject := range subjectNames {
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
							e := esa.MakeEsa(d)
							q := querySeq.Data()
							mat.UpdateMatchLengths(q, e, subjectIDs[j],
								matches.matchLengths[j],
								matches.subjectID[j])
						}
					}
				}
				f.Close()
				matchesSets <- matches
			}(subjectNames, subjectIDs)
		}
		go func() {
			wg.Wait()
			close(matchesSets)
		}()
		matchLengths := make([][]int, 0)
		subjectIDs := make([][]int, 0)
		for _, qs := range querySeqs {
			n := len(qs.Data())
			ml := make([]int, n)
			matchLengths = append(matchLengths, ml)
			n = len(qs.Data())
			ml = make([]int, n)
			subjectIDs = append(subjectIDs, ml)
		}
		for match := range matchesSets {
			for i := 0; i < len(match.matchLengths); i++ {
				for j := 0; j < len(match.matchLengths[i]); j++ {
				}
			}
			for i, lengths := range match.matchLengths {
				for j, length := range lengths {
					println("range lengths: ", j, length)
					if matchLengths[i][j] < length {
						matchLengths[i][j] = length
						subjectIDs[i][j] = match.subjectID[i][j]
					}
				}
			}
		}
		for i, ml := range matchLengths {
			l := 0
			id := 0
			for j := 0; j < len(ml); j++ {
				if ml[j] > l {
					l = ml[j]
					id = subjectIDs[i][j]
				}
				ml[j] = l
				subjectIDs[i][j] = id
				l--
			}
		}
		for i, querySeq := range querySeqs {
			ml := matchLengths[i]
			id := subjectIDs[i]
			f, err = os.Create(querySeq.Header() + ".txt")
			util.Check(err)
			wr := bufio.NewWriter(f)
			var pair []string
			seen := make(map[string]bool)
			for _, ID := range id {
				if seq, exists := sID[ID]; exists {
					if !seen[seq] {
						pair = append(pair,
							fmt.Sprintf("%d=%s", ID, seq))
						seen[seq] = true
					}
				}
			}
			fmt.Fprintf(wr, "#%s\t", query)
			fmt.Fprintf(wr, "%s", strings.Join(pair, "\t"))
			fmt.Fprintf(wr, "\n")
			for j := 0; j < len(ml); j++ {
				fmt.Fprintf(wr, "%d\t%d\n",
					ml[j],
					id[j])
			}
			err = wr.Flush()
			util.Check(err)

			f.Close()
		}
	}
}
