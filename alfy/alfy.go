package main

import (
	"flag"
	"github.com/evolbioinf/alfy/util"
	"github.com/evolbioinf/clio"
)

func main() {
	util.PrepareErrorMessages("alfy")
	optV := flag.Bool("v", false, "version")
	u := "alfy <prepAlfy.out>"
	p := "Calculate homology between queries and subjects"
	e := "alfy alfy.in"
	clio.Usage(u, p, e)
	flag.Parse()
	if *optV {
		util.Version("alfy")
	}

}
