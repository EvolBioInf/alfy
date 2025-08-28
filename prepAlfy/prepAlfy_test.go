package main

import (
	"bytes"
	"os"
	"os/exec"
	"strconv"
	"testing"
)

func TestPrepAlfy(t *testing.T) {
	var tests []*exec.Cmd
	p := "./prepAlfy"
	q := "query"
	s := "subject"
	test := exec.Command(p, "-q", q, "-s", s)
	tests = append(tests, test)
	get, err := test.Output()
	if err != nil {
		t.Error(err)
	}
	f := "r" + strconv.Itoa(1) + ".txt"
	want, err := os.ReadFile(f)
	if err != nil {
		t.Error(err)
	}
	if !bytes.Equal(get, want) {
		t.Errorf("get:\n%s\nwant:\n%s\n", get, want)
	}
}
