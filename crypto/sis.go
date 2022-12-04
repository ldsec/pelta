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

// EmbedSecondary appends a secondary linear relation of form Bs + y = z to this relation.
func (r *LinearRelation) EmbedSecondary(B *fastmath.IntMatrix, y, z *fastmath.IntVec) {
	// Extend the doubles the columns of A with zero to accomodate the new relation's A.
	r.A.ExtendCols(fastmath.NewIntMatrix(r.A.Rows(), y.Size(), r.S.BaseRing()))
	Bp := B.Copy().ExtendCols(fastmath.NewIdIntMatrix(y.Size(), r.S.BaseRing()))
	r.A.ExtendRows(Bp)
	r.S.Append(y)
	r.U.Append(z)
}

func (r *LinearRelation) Copy() LinearRelation {
	return LinearRelation{
		A: r.A.Copy(),
		S: r.S.Copy(),
		U: r.U.Copy(),
	}
}