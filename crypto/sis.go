package crypto

import (
	"github.com/ldsec/codeBase/commitment/fastmath"
)

// SISProblem represents an MSIS problem, i.e., As = u for short s.
type SISProblem struct {
	A *fastmath.IntMatrix
	S *fastmath.IntVec
	U *fastmath.IntVec
}

// NewSISProblem creates a new SIS instance s.t. As = u.
func NewSISProblem(A *fastmath.IntMatrix, s *fastmath.IntVec) SISProblem {
	u := A.MulVec(s)
	return SISProblem{A, s, u}
}

// Rebase splits an SIS problem into multiple SIS problems over a smaller ring.
func (s SISProblem) Rebase(newRing fastmath.RingParams) SISProblem {
	return SISProblem{
		A: s.A.Copy().RebaseRowsLossless(newRing, 0),
		S: s.S.Copy().RebaseLossless(newRing, 0),
		U: s.U.Copy().RebaseLossless(newRing, 0),
	}
}
