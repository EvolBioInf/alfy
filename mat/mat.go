// Package mat implements the matching procedure of alfy.
package mat

import (
	"github.com/evolbioinf/esa"
)

// The function GetMatchLengths takes as arguments a query sequence, the enhanced suffix array (ESA) of a subject sequence, the label of that subject sequence, an array of match lengths, a double array of subject labels, whether or not the query is a reverse strand, and whether or not the program should step through all possible match start positions. It then updates the arrays of match lengths and subject labels.
func GetMatchLengths(q []byte, s *esa.Esa, i int,
	ml []int, su [][]int, rev, step bool) {
	j := 0
	dml := make([]int, len(ml))
	dsu := make([][]int, len(su))
	for i := 0; i < len(ml); i++ {
		dsu[i] = make([]int, 1)
	}
	m := len(q)
	for j < m {
		l := s.MatchPref(q[j:]).L
		p := j
		if rev {
			o := l - 1
			if o < 0 {
				o = 0
			}
			p = m - p - 1 - o
		}
		dml[p] = l
		dsu[p][0] = i
		j++
		if !step {
			j += l
		}
	}
	Update(ml, dml, su, dsu)
	Interpolate(ml, su)
}
func Update(ml, dml []int, su, dsu [][]int) {
	p := 0
	for p < len(ml) {
		if dml[p] > ml[p] {
			ml[p] = dml[p]
			su[p] = su[p][:0]
			su[p] = append(su[p], dsu[p]...)
		} else if dml[p] == ml[p] {
			su[p] = union(su[p], dsu[p])
		}
		p++
	}
}
func union(a, b []int) []int {
	if len(a) > len(b) {
		tmp := a
		a = b
		b = tmp
	}
	dic := make(map[int]bool)
	for _, e := range a {
		dic[e] = true
	}
	for _, e := range b {
		if !dic[e] {
			a = append(a, e)
		}
	}
	return a
}
func Interpolate(ml []int, su [][]int) {
	for i := 1; i < len(ml); i++ {
		prev := ml[i-1] - 1
		curr := ml[i]
		if prev < 0 {
			prev = 0
		}
		if prev > curr {
			ml[i] = prev
			su[i] = su[i][:0]
			su[i] = append(su[i], su[i-1]...)
		} else if prev == curr {
			su[i] = union(su[i], su[i-1])
		}
	}
}
