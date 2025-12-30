package main

import (
	"bufio"
	"flag"
	"fmt"
	"github.com/evolbioinf/alfy/util"
	"github.com/evolbioinf/clio"
	"github.com/evolbioinf/nwk"
	"io"
	"log"
	"slices"
	"sort"
	"strconv"
	"strings"
)

type Tree struct {
	start, end int
	si, ei     int
	neighbors  []string
	mm         int
	query      string
	root       *nwk.Node
}

func parse(r io.Reader, args ...interface{}) {
	query := args[0].(string)
	treeAnalysis := args[1].(bool)
	qi, err := strconv.Atoi(query)
	util.Check(err)
	var nsam, seqLen int
	trees := []*Tree{}
	positions := []float64{}
	haplotypes := []string{}
	start := 0
	sc := bufio.NewScanner(r)
	for sc.Scan() {
		line := sc.Text()
		if len(line) == 0 {
			continue
		}
		fc := line[0]
		lc := line[len(line)-1]
		if fc == 'm' {
			if strings.Index(line, "-T") < 0 {
				log.Fatal("please use tree printing, -T")
			}
			fields := strings.Fields(line)
			nsam, err = strconv.Atoi(fields[1])
			util.Check(err)
			for i := 3; i < len(fields); i++ {
				if fields[i] == "-r" {
					seqLen, err =
						strconv.Atoi(fields[i+2])
					util.Check(err)
				}
			}
		} else if lc == ';' {
			tree := new(Tree)
			arr := strings.Split(line, "]")
			hasExt := false
			if len(arr) == 2 {
				hasExt = true
			}
			ts := arr[0]
			if hasExt {
				ts = arr[1]
			}
			ext := seqLen
			if hasExt {
				var e error
				ext, e = strconv.Atoi(arr[0][1:])
				util.Check(e)
			}
			tree.start = start
			tree.end = tree.start + ext - 1
			tree.query = query
			r := strings.NewReader(ts)
			sc := nwk.NewScanner(r)
			sc.Scan()
			tree.root = sc.Tree()
			start = tree.end + 1
			trees = append(trees, tree)
		} else if fc == 'p' {
			fields := strings.Fields(line)
			for i := 1; i < len(fields); i++ {
				p, e := strconv.ParseFloat(fields[i], 64)
				util.Check(e)
				positions = append(positions, p)
			}
		} else if fc == '1' || fc == '0' {
			haplotypes = append(haplotypes, line)
		}
	}
	if treeAnalysis {
		for _, tree := range trees {
			findParent(tree.root, tree)
		}
	} else {
		si := 0
		for _, tree := range trees {
			np := len(positions)
			s := float64(tree.start) / float64(seqLen)
			for si < np && positions[si] < s {
				si++
			}
			tree.si = si
			ei := si
			e := float64(tree.end) / float64(seqLen)
			for ei < np && positions[ei] <= e {
				ei++
			}
			tree.ei = ei - 1
			si = ei
		}
		dist := make([]int, nsam)
		for _, tree := range trees {
			for i := 0; i < nsam; i++ {
				if i == qi {
					continue
				}
				dist[i] = 0
				for j := tree.si; j <= tree.ei; j++ {
					c1 := haplotypes[i][j]
					c2 := haplotypes[qi][j]
					if c1 != c2 {
						dist[i]++
					}
				}
			}
			min := seqLen + 1
			for i := 0; i < nsam; i++ {
				if i == qi {
					continue
				}
				if dist[i] < min {
					min = dist[i]
				}
			}
			tree.mm = min
			for i := 0; i < nsam; i++ {
				if i+1 != qi && dist[i] == min {
					nei := strconv.Itoa(i + 1)
					tree.neighbors =
						append(tree.neighbors, nei)
				}
			}
		}
	}
	for _, tree := range trees {
		sort.Strings(tree.neighbors)
	}
	c := 0
	for i := 1; i < len(trees); i++ {
		if slices.Equal(trees[c].neighbors,
			trees[i].neighbors) {
			trees[c].end = trees[i].end
			trees[c].mm += trees[i].mm
		} else {
			c++
			trees[c] = trees[i]
		}
	}
	trees = trees[0 : c+1]
	mm := 0.0
	fmt.Printf(">%s\n", query)
	for _, tree := range trees {
		mm = -1
		if !treeAnalysis {
			ext := tree.end - tree.start
			mm = float64(tree.mm) / float64(ext)
		}
		fmt.Printf("%d\t%d\t%.4f\t%s",
			tree.start+1, tree.end+1, mm,
			tree.neighbors[0])
		for _, neighbor := range tree.neighbors[1:] {
			fmt.Printf(" %s", neighbor)
		}
		fmt.Printf("\n")
	}
}
func findParent(v *nwk.Node, tree *Tree) {
	if v != nil {
		if v.Label == tree.query {
			addNeighbors(v.Parent, tree)
		} else {
			findParent(v.Child, tree)
			findParent(v.Sib, tree)
		}
	}
}
func addNeighbors(v *nwk.Node, tree *Tree) {
	if v != nil {
		if v.Child == nil && v.Label != tree.query {
			tree.neighbors = append(tree.neighbors,
				v.Label)
		}
		addNeighbors(v.Child, tree)
		addNeighbors(v.Sib, tree)
	}
}
func main() {
	clio.PrepLog("ms2nn")
	u := "ms2nn -q <query> [option]... [foo.ms]..."
	p := "Convert the output of Hudson's ms to alfy output."
	e := "ms 5 1 -t 100 -r 100 10000 -T | ms2nn -q 1"
	clio.Usage(u, p, e)
	flagV := flag.Bool("v", false, "version")
	flagQ := flag.String("q", "", "query")
	flagT := flag.Bool("t", false, "tree analysis "+
		"(default distance analysis)")
	flag.Parse()
	if *flagV {
		util.Version("ms2nn")
	}
	if *flagQ == "" {
		log.Fatal("please use -q to enter a query")
	}
	files := flag.Args()
	clio.ParseFiles(files, parse, *flagQ, *flagT)
}
