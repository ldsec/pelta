package crypto

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

// ValidityProofDef represents a definition of the proof for the validity of a linear relation As = u.
type ValidityProofDef struct {
	Rel LinearRelation
}

// TernaryProofDef represents a definition of the proof for the ternary structure of a slice over s in a linear relation As = u.
type TernaryProofDef struct {
	Rel    LinearRelation
	Target fastmath.Slice
}

// ApproxBoundProofDef represents a definition of the proof for the approximate infinity bound of a slice over s in a linear relation As = u.
type ApproxBoundProofDef struct {
	Rel    LinearRelation
	Target fastmath.Slice
	Bound  *big.Int
}

// CreateValidityProofDef creates a new proof definition for proving As = u.
func (r LinearRelation) CreateValidityProofDef() ValidityProofDef {
	if !r.IsValid() {
		panic("cannot create a validity proof definition on an invalid linear relation")
	}
	return ValidityProofDef{r}
}

// CreateTernaryProofDef creates a new proof definition for proving the ternary structure of a slice of s.
func (r LinearRelation) CreateTernaryProofDef(start, end int) TernaryProofDef {
	s := fastmath.NewSlice(start, end)
	if !r.S.Slice(s).All(func(el uint64) bool { return el < 3 }) {
		panic("cannot create a ternary structure proof on a non-ternary subvector")
	}
	return TernaryProofDef{Rel: r, Target: s}
}

// CreateABPDef creates a new proof definition for proving an approximate bound on the infinity norm of a slice of s.
func (r LinearRelation) CreateApproxBoundProofDef(start, end int, bound *big.Int) ApproxBoundProofDef {
	s := fastmath.NewSlice(start, end)
	if r.S.Slice(s).Max() > bound.Uint64() {
		panic("cannot create a infinity bound proof on a subvector that has large infinity norm")
	}
	return ApproxBoundProofDef{Rel: r, Target: s, Bound: bound}
}
