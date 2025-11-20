package main

import (
	"bufio"
	"flag"
	"fmt"
	"github.com/evolbioinf/alfy/util"
	"github.com/evolbioinf/clio"
	"log"
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
	start int
	end   int
	score int
	sbp   float64
	ID    []int
}

func main() {
	util.PrepareErrorMessages("alfy")
	optV := flag.Bool("v", false, "version")
	optF := flag.String("f", "", "interval file")
	optW := flag.Int("w", 80, "window length")
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
	sc := bufio.NewScanner(file)
	for sc.Scan() {
		line := sc.Text()
		if len(line) == 0 {
			continue
		}
		if line[0] == '#' {
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
		for _, sequence := range query.sequences {
			m := len(sequence.subjectNames)
			n := sequence.length
			matches := make([][]int, m)
			for i := range sequence.subjectNames {
				matches[i] = make([]int, n)
			}
			for _, interval := range sequence.intervals {
				for _, i := range interval.subjectIDs {
					p := 0
					start := interval.start
					end := interval.start + interval.ml - 1
					for j := start; j <= end; j++ {
						curr := interval.ml - p
						if curr > matches[i][j] {
							matches[i][j] = curr
						}
						p++
					}
				}
			}
			score := make([][]int, len(sequence.subjectNames))
			for i := range matches {
				if len(matches[i]) >= *optW {
					l := 0
					r := 0
					match := 0
					for r < *optW {
						match = match + matches[i][r]
						r++
					}
					score[i] = append(score[i], match)
					match = 0
					for r < len(matches[i]) {
						l++
						r++
						for j := l; j < r; j++ {
							match = match + matches[i][j]
						}
						score[i] = append(score[i], match)
						match = 0
					}
				}
			}
			maxScore := make([]int, 0)
			maxSbjct := make([][]int, 0)
			for j := 0; j < len(score[0]); j++ {
				max := -1
				var ms []int
				for i := 0; i < m; i++ {
					if max < score[i][j] {
						max = score[i][j]
						ms = ms[:0]
						ms = append(ms, i)
					} else if max == score[i][j] {
						ms = append(ms, i)
					}
				}
				maxScore = append(maxScore, max)
				maxSbjct = append(maxSbjct, ms)
				ms = ms[:0]
			}
			windows := []*Window{}
			var window *Window
			p := 0
			s := 0
			for i := 0; i < len(maxSbjct); i++ {
				window = new(Window)
				window.end = *optW + i
				window.ID = maxSbjct[i]
				curr := window.ID
				if i+1 >= len(maxSbjct) {
					window.start = p
					window.score = s + maxScore[i]
					window.sbp = float64(window.score) /
						(float64(window.end) -
							float64(window.start))
					windows = append(windows, window)
					break
				} else {
					next := maxSbjct[i+1]
					if slices.Equal(curr, next) {
						window.end = i + *optW + 1
						s = s + maxScore[i]
					} else {
						window.score = window.score
						window.end = i + *optW
						window.start = p
						window.score = s + maxScore[i]
						window.sbp = float64(window.score) /
							(float64(window.end) - float64(window.start))
						windows = append(windows, window)
						p = i + 1
						s = 0
					}
				}
			}
			fmt.Printf("#Query %s\n>Sequence: %s\n",
				query.name, sequence.name)
			for _, win := range windows {
				str := make([]string, len(win.ID))
				for i, val := range win.ID {
					name := strconv.Itoa(val)
					name = sequence.subjectNames[val]
					e := strings.LastIndex(name, ".")
					name = name[:e]
					str[i] = fmt.Sprintf("%v", name)
				}
				fmt.Printf("%d\t%d\t%d\t%.3f\t%v\n",
					win.start, win.end, win.score,
					win.sbp, strings.Join(str, ","))
			}
		}
	}
}
