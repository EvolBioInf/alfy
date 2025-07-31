package mat

import (
	"github.com/evolbioinf/esa"
	"testing"
)

func check(t *testing.T, get, want []int) {
	for i, _ := range get {
		if get[i] != want[i] {
			t.Errorf("\nget:  %v\n"+
				"want: %v",
				get, want)
		}
	}
}
func TestMat(t *testing.T) {
	queries := [][]byte{
		[]byte("AAATAATATT"),
		[]byte("AAATAGTATT")}
	subjects := [][]byte{
		[]byte("AATTTAAATT"),
		[]byte("TAATAATATT")}
	wml := [][]int{
		[]int{4, 3, 2, 7, 6, 5, 4, 3, 2, 1},
		[]int{4, 3, 2, 2, 1, 0, 4, 3, 2, 1}}
	wsu := [][]int{
		[]int{0, 0, 0, 1, 1, 1, 1, 1, 0, 0},
		[]int{0, 0, 0, 1, 1, 0, 1, 1, 1, 0}}
	for i, q := range queries {
		n := len(q)
		ml := make([]int, n)
		su := make([]int, n)
		for j := 0; j < n; j++ {
			ml[j] = -1
			su[j] = -1
		}
		for j, subject := range subjects {
			s := esa.MakeEsa(subject)
			UpdateMatchLengths(q, s, j, ml, su)
		}
		check(t, ml, wml[i])
		check(t, su, wsu[i])
	}
}
