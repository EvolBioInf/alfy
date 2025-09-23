package mat

import (
	"github.com/evolbioinf/esa"
	"slices"
	"testing"
)

func checkM(t *testing.T, get, want []int) {
	diff := false
	for i, _ := range get {
		if get[i] != want[i] {
			diff = true
			break
		}
	}
	if diff {
		t.Errorf("\nget:  %v\n"+
			"want: %v",
			get, want)
	}
}
func checkS(t *testing.T, get, want [][]int) {
	diff := false
	for i, _ := range get {
		if slices.Compare(get[i], want[i]) != 0 {
			diff = true
			break
		}
	}
	if diff {
		t.Errorf("\nget:  %v\n"+
			"want: %v",
			get, want)
	}
}
func TestMat(t *testing.T) {
	queries := [][]byte{
		[]byte("AAATAATAGT"),
		[]byte("AAATAGTATT")}
	subjects := [][]byte{
		[]byte("AATTTAAATT"),
		[]byte("TAATAATATT")}
	wml := [][]int{
		[]int{4, 3, 2, 5, 4, 3, 2, 1, 0, 1},
		[]int{4, 3, 2, 2, 1, 0, 4, 3, 2, 1}}
	wsu := [][][]int{
		[][]int{{0}, {0}, {0}, {1}, {1},
			{1}, {1}, {1}, {0, 1}, {0, 1}},
		[][]int{{0}, {0}, {0}, {1}, {1},
			{0, 1}, {1}, {1}, {1}, {0, 1}}}
	rev := false
	for i, q := range queries {
		n := len(q)
		ml := make([]int, n)
		su := make([][]int, n)
		for j := 0; j < n; j++ {
			ml[j] = -1
			su[j] = append(su[j], -1)
		}
		for j, subject := range subjects {
			s := esa.MakeEsa(subject)
			GetMatchLengths(q, s, j, ml, su, rev)
		}
		checkM(t, ml, wml[i])
		checkS(t, su, wsu[i])
	}
	queries = [][]byte{
		[]byte("ACTATTATTT"),
		[]byte("AATACTATTT")}
	wml = [][]int{
		[]int{2, 2, 1, 2, 4, 3, 2, 1, 0, 1},
		[]int{2, 4, 3, 2, 1, 0, 4, 3, 2, 1}}
	wsu = [][][]int{
		[][]int{{0}, {1}, {1}, {0}, {1},
			{1}, {0, 1}, {0, 1}, {0, 1}, {0, 1}},
		[][]int{{0}, {1}, {1}, {0, 1}, {0, 1},
			{0, 1}, {1}, {0, 1}, {0, 1}, {0, 1}}}
	rev = true
	for i, q := range queries {
		n := len(q)
		ml := make([]int, n)
		su := make([][]int, n)
		for j := 0; j < n; j++ {
			ml[j] = -1
			su[j] = append(su[j], -1)
		}
		for j, subject := range subjects {
			s := esa.MakeEsa(subject)
			GetMatchLengths(q, s, j, ml, su, rev)
		}
		checkM(t, ml, wml[i])
		checkS(t, su, wsu[i])
	}
}
