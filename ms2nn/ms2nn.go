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
	neighbors  []string
	query      string
	root       *nwk.Node
}

func parse(r io.Reader, args ...interface{}) {
	query := args[0].(string)
	mergeT := args[1].(bool)
	mergeD := args[2].(bool)
	trees := []*Tree{}
	sc := bufio.NewScanner(r)
	start := 1
	for sc.Scan() {
		line := sc.Text()
		if len(line) > 0 && line[0] == '[' {
			tree := new(Tree)
			arr := strings.Split(line, "]")
			ts := arr[1]
			ext, e := strconv.Atoi(arr[0][1:])
			util.Check(e)
			tree.start = start
			tree.end = tree.start + ext - 1
			tree.query = query
			r := strings.NewReader(ts)
			sc := nwk.NewScanner(r)
			sc.Scan()
			tree.root = sc.Tree()
			start = tree.end + 1
			trees = append(trees, tree)
		} else {
			continue
		}
	}
	for _, tree := range trees {
		findParent(tree.root, tree)
	}
	for _, tree := range trees {
		sort.Strings(tree.neighbors)
	}
	if mergeT {
		c := 0
		for i := 1; i < len(trees); i++ {
			if slices.Equal(trees[c].neighbors, trees[i].neighbors) {
				trees[c].end = trees[i].end
			} else {
				c++
				trees[c] = trees[i]
			}
		}
		trees = trees[0 : c+1]

	} else if mergeD {

	}
	fmt.Printf(">%s\n", query)
	for _, tree := range trees {
		fmt.Printf("%d\t%d\t-1\t%s",
			tree.start, tree.end, tree.neighbors[0])
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
			tree.neighbors = append(tree.neighbors, v.Label)
		}
		addNeighbors(v.Child, tree)
		addNeighbors(v.Sib, tree)
	}
}
func main() {
	clio.PrepLog("ms2nn")
	u := "ms2nn [-h -v] -q <query> [foo.ms]..."
	p := "Convert the output of Hudson's ms to alfy output."
	e := "ms 5 1 -t 100 -r 100 10000 -T | ms2nn -q 3"
	clio.Usage(u, p, e)
	flagV := flag.Bool("v", false, "version")
	flagQ := flag.String("q", "", "query")
	flagT := flag.Bool("t", false, "merge by topology")
	flagD := flag.Bool("d", false, "merge by distance")
	flag.Parse()
	if *flagV {
		util.Version("ms2nn")
	}
	if *flagQ == "" {
		log.Fatal("please use -q to enter a query")
	}
	if *flagT && *flagD {
		m := "please use at most one of the two " +
			"merge flags -t and -d"
		log.Fatal(m)
	}
	files := flag.Args()
	clio.ParseFiles(files, parse, *flagQ, *flagT, *flagD)
}
