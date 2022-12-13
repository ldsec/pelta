package crypto

import (
	"fmt"

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

func NewLinearRelationWithLHS(A *fastmath.IntMatrix, s, u *fastmath.IntVec) LinearRelation {
	return LinearRelation{A, s, u}
}

func (r LinearRelation) Rebase(newRing fastmath.RingParams) LinearRelation {
	return LinearRelation{
		A: r.A.Copy().RebaseRowsLossless(newRing, 0),
		S: r.S.Copy().RebaseLossless(newRing, 0),
		U: r.U.Copy().RebaseLossless(newRing, 0),
	}
}

// AppendDependent appends a relation dependent on s, i.e., B(s, y) = (u, z), by performing necessary zero-padding on A.
// The resulting relation is (A || 0, B) (s, y) = (u, z)
func (r *LinearRelation) AppendDependent(B *fastmath.IntMatrix, y, z *fastmath.IntVec) *LinearRelation {
	n := r.A.Cols()
	m := r.A.Rows()
	np := B.Cols()
	if n > np {
		panic("cannot append dependent relation with cols(A) > cols(B)")
	} else if np > n {
		r.A.ExtendCols(fastmath.NewIntMatrix(m, np-n, r.S.BaseRing()))
	}
	r.A.ExtendRows(B)
	if np > n {
		r.S.Append(y)
	}
	r.U.Append(z)
	return r
}

// AppendIndependent appends an independent linear relation of form By = z to this relation As = u.
// The resulting relation is (A || 0, 0 || B) (s, y) = (u, z).
func (r *LinearRelation) AppendIndependent(rp LinearRelation) *LinearRelation {
	m := r.A.Rows()
	n := r.A.Cols()
	mp := rp.A.Rows()
	np := rp.A.Cols()
	AHorizontalExt := fastmath.NewIntMatrix(m, np, r.S.BaseRing())
	r.A.ExtendCols(AHorizontalExt)
	AVerticalExt := fastmath.NewIntMatrix(mp, n, r.S.BaseRing()).ExtendCols(rp.A)
	r.A.ExtendRows(AVerticalExt)
	r.S.Append(rp.S)
	r.U.Append(rp.U)
	return r
}

// Extend extends a linear relation As = u with By = z s.t. As + By = z + u.
// The resulting relation is (A || B) (s, y) = z + u
func (r *LinearRelation) Extend(rp LinearRelation) *LinearRelation {
	r.ExtendPartial(rp.A, rp.S)
	r.U.Add(rp.U)
	return r
}

// ExtendPartial extends a linear relation As = u with (B, y) s.t. As + By = u.
// The resulting relation is (A || B) (s, y) = z
func (r *LinearRelation) ExtendPartial(B *fastmath.IntMatrix, y *fastmath.IntVec) *LinearRelation {
	m := r.A.Rows()
	mp := B.Rows()
	if mp != m {
		panic("cannot extend (B, y) because B has incompatible number of rows")
	}
	r.A.ExtendCols(B)
	r.S.Append(y)
	return r
}

// AppendDependentOnS appends a secondary linear relation of form Bs + y = z to this relation As = u.
// Note: B must be a square matrix with rows(B) = cols(B) = cols(A)
func (r *LinearRelation) AppendDependentOnS(B *fastmath.IntMatrix, y, z *fastmath.IntVec) *LinearRelation {
	n := r.A.Cols()
	mp := B.Rows()
	np := B.Cols()
	if np != n {
		panic("cannot append Bs + y = z where B has different column size")
	}
	BExtended := B.Copy().ExtendCols(fastmath.NewIdIntMatrix(mp, r.S.BaseRing()))
	return r.AppendDependent(BExtended, y, z)
}

func (r *LinearRelation) Copy() LinearRelation {
	return LinearRelation{
		A: r.A.Copy(),
		S: r.S.Copy(),
		U: r.U.Copy(),
	}
}

func (r *LinearRelation) String() string {
	return fmt.Sprintf("A: %s\ns: %s\nu: %s", r.A.String(), r.S.String(), r.U.String())
}

func (r *LinearRelation) SizesString() string {
	return fmt.Sprintf("A[%dx%d] s[%d] = u[%d]", r.A.Rows(), r.A.Cols(), r.S.Size(), r.U.Size())
}

// Verify returns true iff As = u holds.
func (r *LinearRelation) Verify() bool {
	return r.A.MulVec(r.S).Eq(r.U)
}
