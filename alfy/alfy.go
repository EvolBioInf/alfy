package main

import (
	"bufio"
	"flag"
	"fmt"
	"github.com/evolbioinf/alfy/util"
	"github.com/evolbioinf/clio"
	"github.com/evolbioinf/sus"
	"log"
	"math"
	"os"
	"slices"
	"strconv"
	"strings"
)

type Query struct {
	name      string
	sequences []*Sequence
}
type Sequence struct {
	name         string
	length       int
	subjectNames map[int]string
	intervals    []*Interval
}
type Interval struct {
	start      int
	end        int
	ml         int
	subjectIDs []int
}
type Window struct {
	start  int
	end    int
	score  []int
	winner []int
}
type Summary struct {
	start int
	end   int
	nm    float64
	label []int
}

func main() {
	util.PrepareErrorMessages("alfy")
	optV := flag.Bool("v", false, "version")
	optF := flag.String("f", "", "interval file")
	optW := flag.Int("w", 80, "window length")
	optQ := flag.Float64("q", 0.1, "quantile of match length distribution")
	u := "alfy [prepAlfy.out]"
	p := "Calculate homology between queries and subjects."
	e := "alfy alfy.in"
	clio.Usage(u, p, e)
	flag.Parse()
	if *optV {
		util.Version("alfy")
	}
	if *optF == "" {
		m := "please provide the intervals file " +
			"(prepAlfy output)"
		fmt.Fprintf(os.Stderr, "%s\n", m)
		os.Exit(1)
	}
	file, err := os.Open(*optF)
	util.Check(err)
	defer file.Close()
	queries := []*Query{}
	var query *Query
	var sequence *Sequence
	var interval *Interval
	var slen int
	var gc float64
	sc := bufio.NewScanner(file)
	for sc.Scan() {
		line := sc.Text()
		if len(line) == 0 {
			continue
		}
		if line[0] == '-' {
			var err error
			fields := strings.Fields(line)
			slen, err = strconv.Atoi(fields[1])
			if err != nil {
				log.Fatalf("%q is not a number", fields[1])
			}
			gc, err = strconv.ParseFloat(fields[2], 64)
			if err != nil {
				log.Fatalf("%q is not a float", fields[2])
			}
		} else if line[0] == '#' {
			query = new(Query)
			query.name = line[1:]
			queries = append(queries, query)
		} else if line[0] == '>' {
			var err error
			sequence = new(Sequence)
			fields := strings.Fields(line)
			sequence.name = fields[0][1:]
			sequence.length, err = strconv.Atoi(fields[1])
			if err != nil {
				log.Fatalf("%q is not a number", fields[1])
			}
			if sequence.subjectNames == nil {
				sequence.subjectNames = make(map[int]string)
			}
			for i := 2; i < len(fields); i++ {
				arr := strings.Split(fields[i], "=")
				k, err := strconv.Atoi(arr[0])
				if err != nil {
					log.Fatalf("%q is not a number",
						arr[0])
				}
				sequence.subjectNames[k-1] = arr[1]
			}
			query.sequences = append(query.sequences, sequence)
		} else {
			var err error
			interval = new(Interval)
			fields := strings.Fields(line)
			interval.start, err = strconv.Atoi(fields[0])
			if err != nil {
				log.Fatalf("%q is not a number",
					fields[0])
			}
			interval.end, err = strconv.Atoi(fields[1])
			if err != nil {
				log.Fatalf("%q is not a number",
					fields[1])
			}
			interval.ml, err = strconv.Atoi(fields[2])
			if err != nil {
				log.Fatalf("%q is not a number",
					fields[2])
			}
			arr := strings.Split(fields[3], ",")
			for i := 0; i < len(arr); i++ {
				x, err := strconv.Atoi(arr[i])
				if err != nil {
					log.Fatalf("%q is not a number ", arr[i])
				}
				interval.subjectIDs =
					append(interval.subjectIDs, x-1)
			}
			sequence.intervals = append(sequence.intervals,
				interval)
		}
	}
	for _, query := range queries {
		for i, sequence := range query.sequences {
			m := len(sequence.subjectNames)
			n := sequence.length
			matches := make([][]int, m)
			for a := range sequence.subjectNames {
				matches[a] = make([]int, n)
			}
			for _, interval := range sequence.intervals {
				for _, a := range interval.subjectIDs {
					p := 0
					start := interval.start
					end := interval.start + interval.ml - 1
					for b := start; b <= end; b++ {
						curr := interval.ml - p
						if curr > matches[a][b] {
							matches[a][b] = curr
						}
						p++
					}
				}
			}
			max := make([]int, len(matches[0]))
			maxID := make([][]int, len(matches[0]))
			for b := range matches[0] {
				for a := range matches {
					if matches[a][b] > max[b] {
						max[b] = matches[a][b]
						maxID[b] = maxID[b][:0]
						maxID[b] = append(maxID[b], a)
					} else if matches[a][b] == max[b] {
						maxID[b] = append(maxID[b], a)
					}
				}
			}
			mismatches := make([]int, len(matches[0]))
			for a := 0; a < len(max); a++ {
				if a == 0 {
					continue
				} else if max[a] >= max[a-1] && max[a] != 0 {
					mismatches[a] = 1
				}
			}

			q := sus.Quantile(slen, gc, *optQ) - 1
			t := int(math.Round(float64(*optW) / float64(q)))
			score := make([]int, len(sequence.subjectNames))
			windows := make([][]*Window, len(query.sequences))
			var window *Window
			if len(max) >= *optW {
				nm := 0
				l := 0
				r := 0
				for r < *optW {
					if mismatches[r] == 1 {
						nm++
					}
					for _, Id := range maxID[r] {
						score[Id]++
					}
					r++
				}
				open := false
				prev := maxID[r-1]
				shifted := false
				for r < len(max) {
					shifted = true
					if nm < t {
						if open {
							if slices.Equal(prev, maxID[r]) {
								window.end = r
							} else {
								windows[i] = append(windows[i], window)
								window = new(Window)
								window.score = make([]int, len(score))
								copy(window.score, score)
								window.start = l
								window.end = r
							}
						} else {
							window = new(Window)
							window.score = make([]int, len(score))
							copy(window.score, score)
							window.start = l
							window.end = r
							open = true
						}
					} else if open && window.end < l {
						open = false
						windows[i] = append(windows[i], window)
					}
					if mismatches[l] == 1 {
						nm--
					}
					if mismatches[r] == 1 {
						nm++
					}
					for _, a := range maxID[l] {
						score[a]--
					}
					if r < len(max)-1 {
						for _, a := range maxID[r] {
							score[a]++
						}
					}
					prev = maxID[r]
					l++
					r++
				}
				if !shifted && r == len(max) {
					if nm < t {
						window = new(Window)
						window.score = make([]int, len(score))
						copy(window.score, score)
						window.start = l
						window.end = r - 1
						windows[i] = append(windows[i], window)
					}
				}
				if open {
					windows[i] = append(windows[i], window)
				}
			}
			for _, win := range windows {
				for _, w := range win {
					max := 0
					sID := make([]int, 0)
					for a, count := range w.score {
						if count > max {
							max = count
							sID = sID[:0]
							sID = append(sID, a)
						} else if count == max {
							sID = append(sID, a)
						}
					}
					w.winner = make([]int, len(sID))
					copy(w.winner, sID)
				}
			}
			var summary []*Summary
			var prev *Summary
			for _, win := range windows {
				for _, w := range win {
					if prev == nil {
						prev = &Summary{
							start: w.start,
							end:   w.end,
							label: w.winner}
					} else if slices.Equal(w.winner,
						prev.label) {
						prev.end = w.end
					} else {
						prev = &Summary{
							start: w.start,
							end:   w.end,
							label: w.winner}
						summary = append(summary,
							prev)
					}
				}
			}
			if prev != nil {
				summary = append(summary, prev)
			}
			for _, sum := range summary {
				var len int
				var nm int
				len = sum.end - sum.start + 1
				for a := sum.start; a <= sum.end; a++ {
					if mismatches[a] == 1 {
						nm++
					}
				}
				sum.nm = float64(nm) / float64(len)
			}
			fmt.Printf("#Query %s\n>Sequence: %s\n",
				query.name, sequence.name)
			for _, sum := range summary {
				str := make([]string, len(sum.label))
				for i, val := range sum.label {
					var name string
					if val < 0 {
						name = "non-homologous"
					} else {
						name = strconv.Itoa(val)
						name = sequence.subjectNames[val]
						e := strings.LastIndex(name, ".")
						name = name[:e]
					}
					str[i] = fmt.Sprintf("%v", name)
				}
				fmt.Printf("%d\t%d\t%.3f\t%v\n",
					sum.start+1, sum.end+1,
					sum.nm, strings.Join(str, ","))
			}
		}
	}
}
