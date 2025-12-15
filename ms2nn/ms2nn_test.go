package main

import (
	"bytes"
	"os"
	"os/exec"
	"strconv"
	"testing"
)

func TestMs2nn(t *testing.T) {
	var tests []*exec.Cmd
	p := "./ms2nn"
	f := "msTest.out"
	test := exec.Command(p, "-q", "1", f)
	tests = append(tests, test)
	test = exec.Command(p, "-q", "2", f)
	tests = append(tests, test)
	test = exec.Command(p, "-q", "1", "-t", f)
	tests = append(tests, test)
	test = exec.Command(p, "-q", "2", "-t", f)
	tests = append(tests, test)
	for i, test := range tests {
		get, err := test.Output()
		if err != nil {
			t.Error(err)
		}
		f = "r" + strconv.Itoa(i+1) + ".txt"
		want, err := os.ReadFile(f)
		if err != nil {
			t.Error(err)
		}
		if !bytes.Equal(get, want) {
			t.Errorf("get:\n%s\nwant:\n%s\n", get, want)
		}
	}
}
