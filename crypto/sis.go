package crypto

import (
	"github.com/ldsec/codeBase/commitment/fastmath"
)

// LinearRelation represents a linear relation, i.e., As = u for some s.
type LinearRelation struct {
	A *fastmath.IntMatrix
	S *fastmath.IntVec
	U *fastmath.IntVec
}

// NewLinearRelation creates a new linear relation instance s.t. As = u.
func NewLinearRelation(A *fastmath.IntMatrix, s *fastmath.IntVec) LinearRelation {
	u := A.MulVec(s)
	return LinearRelation{A, s, u}
}

func (r LinearRelation) Rebase(newRing fastmath.RingParams) LinearRelation {
	return LinearRelation{
		A: r.A.Copy().RebaseRowsLossless(newRing, 0),
		S: r.S.Copy().RebaseLossless(newRing, 0),
		U: r.U.Copy().RebaseLossless(newRing, 0),
	}
}
