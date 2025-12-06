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
	trees := []*Tree{}
	prev := ""
	var tree *Tree
	start := 1
	sc := bufio.NewScanner(r)
	for sc.Scan() {
		line := sc.Text()
		if len(line) > 0 && line[0] == '[' {
			arr := strings.Split(line, "]")
			curr := arr[1]
			ext, e := strconv.Atoi(arr[0][1:])
			util.Check(e)
			if curr != prev {
				tree = new(Tree)
				tree.start = start
				tree.end = tree.start + ext - 1
				tree.query = query
				r := strings.NewReader(curr)
				sc := nwk.NewScanner(r)
				sc.Scan()
				tree.root = sc.Tree()
				trees = append(trees, tree)
				start = tree.end + 1
			} else {
				tree.end += ext
				start = tree.end + 1
			}
		} else {
			continue
		}
	}
	for _, tree := range trees {
		findParent(tree.root, tree)
	}
	for _, tree := range trees {
		if len(tree.neighbors) == 0 {
			fmt.Println(tree)
		}
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
	clio.PrepLog("ms2alfy")
	u := "ms2alfy [-h -v] -q <query> [foo.ms]..."
	p := "Convert the output of Hudson's ms to alfy output."
	e := "ms 5 1 -t 100 -r 100 10000 -T | ms2alfy -q 3"
	clio.Usage(u, p, e)
	flagV := flag.Bool("v", false, "version")
	flagQ := flag.String("q", "", "query")
	flag.Parse()
	if *flagV {
		util.Version("ms2alfy")
	}
	if *flagQ == "" {
		log.Fatal("please use -q to enter a query")
	}
	files := flag.Args()
	clio.ParseFiles(files, parse, *flagQ)
}
