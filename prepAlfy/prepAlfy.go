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
	optT := flag.Int("t", ncpu, "number of threads")
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
		for i, Esa := range EsaSets {
			subjectIDs := subjectIDsets[i]
			wg.Add(1)
			go func(Esa []*esa.Esa, subjectIDs []int) {
				defer wg.Done()
				var matches Matches
				n := len(querySeqs)
				matches.matchLengths = make([][]int, n)
				matches.subjectID = make([][]int, n)
				for i, querySeq := range querySeqs {
					m := len(querySeq.Data())
					matches.matchLengths[i] = make([]int, m)
					matches.subjectID[i] = make([]int, m)
					for j, esa := range Esa {
						q := querySeq.Data()
						//  println(esa)
						mat.UpdateMatchLengths(q, esa, subjectIDs[j],
							matches.matchLengths[i],
							matches.subjectID[i])
					}
				}
				matchesSets <- matches
			}(Esa, subjectIDs)
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
					if matchLengths[i][j] < length {
						matchLengths[i][j] = length
						subjectIDs[i][j] =
							match.subjectID[i][j]
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
		for i := range querySeqs {
			ml := matchLengths[i]
			id := subjectIDs[i]
			wr := bufio.NewWriter(os.Stdout)
			var pair []string
			seen := make(map[string]bool)
			for _, ID := range id {
				if seq, exists := sID[ID]; exists {
					if !seen[seq] {
						pair = append(pair,
							fmt.Sprintf("%d=%s",
								ID+1, seq))
						seen[seq] = true
					}
				}
			}

			sort.Strings(pair)
			fmt.Fprintf(wr, "#%s\t", query)
			fmt.Fprintf(wr, "%s", strings.Join(pair, "\t"))
			fmt.Fprintf(wr, "\n")
			for j := 0; j < len(ml); j++ {
				fmt.Fprintf(wr, "%d\t%d\n",
					ml[j],
					id[j]+1)
			}
			wr.Flush()
		}
	}
}
