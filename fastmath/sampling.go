package fastmath

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// PolySampler represents a random polynomial sampler.
type PolySampler interface {
	Read(pol *ring.Poly)
}

// NewRandomPolynomial returns a random polynomial sampled from the given `sampler`.
func NewRandomPolynomial(sampler PolySampler, baseRing *ring.Ring) Poly {
	g := NewZeroPoly(baseRing)
	sampler.Read(g.ref)
	return g
}

// NewRandomPolynomialVector constructs a vector, whose elements sampled from the given `sampler`.
func NewRandomPolynomialVector(size int, sampler PolySampler, baseRing *ring.Ring) PolyVec {
	v := NewPolyVec(size, baseRing)
	v.Populate(func(i int) Poly {
		return NewRandomPolynomial(sampler, baseRing)
	})
	return v
}

// NewRandomPolynomialMatrix constructs a 2D matrix, whose elements sampled from the given `sampler`.
func NewRandomPolynomialMatrix(rows int, cols int, sampler PolySampler, baseRing *ring.Ring) PolyMatrix {
	A := NewPolyMatrix(rows, cols, baseRing)
	A.PopulateRows(func(_ int) PolyVec {
		return NewRandomPolynomialVector(cols, sampler, baseRing)
	})
	return A
}

// NewRandomIntegerVector constructs a random vector of integers mod n.
func NewRandomIntegerVector(size int, n *big.Int, baseRing *ring.Ring) IntVec {
	v := NewIntVec(size, baseRing)
	v.Populate(func(_ int) uint64 {
		return ring.RandInt(n).Uint64()
	})
	return v
}

// NewRandomTernaryIntegerVector constructs a random vector of integers where each element \in {0, 1, 2}.
func NewRandomTernaryIntegerVector(size int, baseRing *ring.Ring) IntVec {
	return NewRandomIntegerVector(size, big.NewInt(3), baseRing)
}

// NewRandomIntegerMatrix constructs a random 2D matrix of integers mod n.
func NewRandomIntegerMatrix(rows int, cols int, n *big.Int, baseRing *ring.Ring) IntMatrix {
	A := NewIntMatrix(rows, cols, baseRing)
	A.PopulateRows(func(_ int) IntVec {
		return NewRandomIntegerVector(cols, n, baseRing)
	})
	return A
}
