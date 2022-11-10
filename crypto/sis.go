package crypto

import "github.com/ldsec/codeBase/commitment/fastmath"

// SISProblem represents an MSIS problem, i.E., As = U for short S.
type SISProblem struct {
	A *fastmath.IntMatrix
	S *fastmath.IntVec
	U *fastmath.IntVec
}

// NewSISProblem creates a new SIS instance S.t. As = u.
func NewSISProblem(A *fastmath.IntMatrix, s *fastmath.IntVec) SISProblem {
	u := A.MulVec(s)
	return SISProblem{A, s, u}
}
