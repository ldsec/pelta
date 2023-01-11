package fastmath

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// PolySampler represents a random polynomial sampler.
type PolySampler interface {
	Read(pol *ring.Poly)
}

// AugmentedTernarySampler samples from {0, 1, 2} instead of {-1, 0, 1}
type AugmentedTernarySampler struct {
	ternarySampler PolySampler
	baseRing       *ring.Ring
}

func (s *AugmentedTernarySampler) Read(pol *ring.Poly) {
	// Read -1, 0, 1
	s.ternarySampler.Read(pol)
	// Convert into 0, 1, 2
	s.baseRing.AddScalar(pol, 1, pol)
}

func NewAugmentedTernarySampler(ternarySampler PolySampler, baseRing *ring.Ring) PolySampler {
	return &AugmentedTernarySampler{ternarySampler, baseRing}
}

// NewRandomPoly returns a random polynomial sampled from the given `sampler`.
func NewRandomPoly(sampler PolySampler, baseRing *ring.Ring) *Poly {
	g := NewPoly(baseRing)
	sampler.Read(g.ref)
	g.SetDirty()
	return g
}

func NewRandomTernaryPoly(baseRing *ring.Ring) *Poly {
	g := NewPoly(baseRing)
	for i := 0; i < g.N(); i++ {
		rand := ring.RandInt(big.NewInt(3)).Uint64()
		g.Set(i, rand)
	}
	g.SetDirty()
	return g
}

// NewRandomPolyVec constructs a vector, whose elements sampled from the given `sampler`.
func NewRandomPolyVec(size int, sampler PolySampler, baseRing *ring.Ring) *PolyVec {
	v := NewPolyVec(size, baseRing)
	v.Populate(func(i int) *Poly {
		return NewRandomPoly(sampler, baseRing)
	})
	return v
}

// NewRandomTernaryPolyVec constructs a vector of ternary polynomials.
func NewRandomTernaryPolyVec(size int, baseRing *ring.Ring) *PolyVec {
	v := NewPolyVec(size, baseRing)
	v.Populate(func(i int) *Poly {
		return NewRandomTernaryPoly(baseRing)
	})
	return v
}

// NewRandomTernaryPolyVec constructs a vector of ternary polynomials.
func NewRandomTernaryPolyMatrix(rows, cols int, baseRing *ring.Ring) *PolyMatrix {
	A := NewPolyMatrix(rows, cols, baseRing)
	A.PopulateRows(func(_ int) *PolyVec {
		return NewRandomTernaryPolyVec(cols, baseRing)
	})
	return A
}

// NewRandomPolyMatrix constructs a 2D matrix, whose elements sampled from the given `sampler`.
func NewRandomPolyMatrix(rows int, cols int, sampler PolySampler, baseRing *ring.Ring) *PolyMatrix {
	A := NewPolyMatrix(rows, cols, baseRing)
	A.PopulateRows(func(_ int) *PolyVec {
		return NewRandomPolyVec(cols, sampler, baseRing)
	})
	return A
}

// NewRandomIntVecFast constructs a random vector of integers from the given poly sampler.
func NewRandomIntVecFast(size int, sampler PolySampler, baseRing *ring.Ring) *IntVec {
	v := NewIntVec(size, baseRing)
	numPolys := len(v.UnderlyingPolys())
	randomPolys := make([]*Poly, numPolys)
	for i := 0; i < numPolys; i++ {
		randomPolys[i] = NewRandomPoly(sampler, baseRing)
	}
	v.SetUnderlyingPolys(randomPolys)
	return v
}

// NewRandomIntVec constructs a random vector of integers mod n.
func NewRandomIntVec(size int, n *big.Int, baseRing *ring.Ring) *IntVec {
	v := NewIntVec(size, baseRing)
	v.Populate(func(_ int) uint64 {
		return ring.RandInt(n).Uint64()
	})
	return v
}

// NewRandomTernaryIntVec constructs a random vector of integers where each element \in {0, 1, 2}.
func NewRandomTernaryIntVec(size int, baseRing *ring.Ring) *IntVec {
	return NewRandomIntVec(size, big.NewInt(3), baseRing)
}

// NewRandomTernaryIntVec constructs a random vector of integers where each element \in {0, 1, 2}.
func NewRandomBinaryIntVec(size int, baseRing *ring.Ring) *IntVec {
	return NewRandomIntVec(size, big.NewInt(2), baseRing)
}

// NewRandomIntMatrixFast constructs a random 2D matrix of integers from the given sampler.
func NewRandomIntMatrixFast(rows int, cols int, sampler PolySampler, baseRing *ring.Ring) *IntMatrix {
	A := NewIntMatrix(rows, cols, baseRing)
	A.PopulateRows(func(_ int) *IntVec {
		return NewRandomIntVecFast(cols, sampler, baseRing)
	})
	return A
}

// NewRandomIntMatrix constructs a random 2D matrix of integers mod n.
func NewRandomIntMatrix(rows int, cols int, n *big.Int, baseRing *ring.Ring) *IntMatrix {
	A := NewIntMatrix(rows, cols, baseRing)
	A.PopulateRows(func(_ int) *IntVec {
		return NewRandomIntVec(cols, n, baseRing)
	})
	return A
}

// NewRandomTernaryIntMatrix constructs a random 2D matrix of integers mod 3.
func NewRandomTernaryIntMatrix(rows int, cols int, baseRing *ring.Ring) *IntMatrix {
	A := NewIntMatrix(rows, cols, baseRing)
	A.PopulateRows(func(_ int) *IntVec {
		return NewRandomTernaryIntVec(cols, baseRing)
	})
	return A
}

// NewRandomBinaryIntMatrix constructs a random 2D matrix of integers mod 2.
func NewRandomBinaryIntMatrix(rows int, cols int, baseRing *ring.Ring) *IntMatrix {
	A := NewIntMatrix(rows, cols, baseRing)
	A.PopulateRows(func(_ int) *IntVec {
		return NewRandomBinaryIntVec(cols, baseRing)
	})
	return A
}
