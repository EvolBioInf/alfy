// Package mat implements the matching procedure of alfy.
package mat

import (
	"github.com/evolbioinf/esa"
)

// The function UpdateMatchLengths takes as arguments a query sequence, the enhanced suffix array (ESA) of a subject sequence, the label of that subject sequence, an array of match lengths, an array of subject labels, and whether or not the query is a reverse strand. It then updates the arrays of match lengths and subject labels.
func UpdateMatchLengths(q []byte, s *esa.Esa, i int,
	ml, su []int, rev bool) {
	j := 0
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
		if l > ml[p] {
			ml[p] = l
			su[p] = i
		}
		j += l + 1
	}
	for i := 1; i < len(q); i++ {
		x := ml[i-1] - 1
		if x >= ml[i] {
			ml[i] = x
			su[i] = su[i-1]
		}
	}
}
