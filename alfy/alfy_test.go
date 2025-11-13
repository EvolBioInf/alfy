package main

import (
	"bytes"
	"os"
	"os/exec"
	"strconv"
	"testing"
)

func TestAlfy(t *testing.T) {
	var tests []*exec.Cmd
	p := "./alfy"
	f := "prep.out"
	w := "5"
	test := exec.Command(p, "-f", f, "-w", w)
	tests = append(tests, test)
	w = "10"
	test = exec.Command(p, "-f", f, "-w", w)
	tests = append(tests, test)
	for i, test := range tests {
		get, err := test.Output()
		if err != nil {
			t.Error(err)
		}
		r := "r" + strconv.Itoa(i+1) + ".txt"
		want, err := os.ReadFile(r)
		if err != nil {
			t.Error(err)
		}
		if !bytes.Equal(get, want) {
			t.Errorf("get:\n%s\nwant:\n%s\n", get, want)
		}
	}
}
