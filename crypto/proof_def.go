package crypto

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

// LinearRelationProofDef represents a definition of the proof for the validity of a linear relation As = u.
type LinearRelationProofDef struct {
	rel LinearRelation
}

// TernaryProofDef represents a definition of the proof for the ternary structure of a slice over some vector.
type TernaryProofDef struct {
	target    fastmath.Slice
	targetVec *fastmath.IntVec
}

// ABPDef represents a definition of the proof for the approximate infinity bound of a slice over some vector.
type ABPDef struct {
	target    fastmath.Slice
	targetVec *fastmath.IntVec
	bound     *big.Int
}

// CreateValidityProofDef creates a new proof definition for proving As = u.
func (r LinearRelation) CreateValidityProofDef() LinearRelationProofDef {
	if !r.IsValid() {
		panic("cannot create a validity proof definition on an invalid linear relation")
	}
	return LinearRelationProofDef{r}
}

// CreateTernaryProofDef creates a new proof definition for proving the ternary structure of a slice of s.
func (r LinearRelation) CreateTernaryProofDef(start, end int) TernaryProofDef {
	s := fastmath.NewSlice(start, end)
	if !r.S.Slice(s).All(func(el uint64) bool { return el < 3 }) {
		panic("cannot create a ternary structure proof on a non-ternary subvector")
	}
	return TernaryProofDef{target: s, targetVec: r.S}
}

// CreateABPDef creates a new proof definition for proving an approximate bound on the infinity norm of a slice of u.
func (r LinearRelation) CreateABPDef(start, end int, bound *big.Int) ABPDef {
	s := fastmath.NewSlice(start, end)
	if r.U.Slice(s).Max() > bound.Uint64() {
		panic("cannot create a infinity bound proof on a subvector that has large infinity norm")
	}
	return ABPDef{target: s, targetVec: r.U, bound: bound}
}
