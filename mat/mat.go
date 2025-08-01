// Package mat implements the matching procedure of alfy.
package mat

import (
	"github.com/evolbioinf/esa"
)

// The function UpdateMatchLengths takes as arguments a query sequence, the enhanced suffix array (ESA) of a subject sequence, the label of that subject sequence, an array of match lengths, and an array of subject labels. It then updates the arrays of match lengths and subject labels.
func UpdateMatchLengths(q []byte, s *esa.Esa, i int,
	ml, su []int) {
	j := 0
	for j < len(q) {
		l := s.MatchPref(q[j:]).L
		if l > ml[j] {
			ml[j] = l
			su[j] = i
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
