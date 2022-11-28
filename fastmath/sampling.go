package fastmath

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// PolySampler represents a random polynomial sampler.
type PolySampler interface {
	Read(pol *ring.Poly)
}

// NewRandomPoly returns a random polynomial sampled from the given `sampler`.
func NewRandomPoly(sampler PolySampler, baseRing *ring.Ring) *Poly {
	g := NewZeroPoly(baseRing)
	sampler.Read(g.ref)
	return g
}

func NewRandomTernaryPoly(baseRing *ring.Ring) *Poly {
	g := NewZeroPoly(baseRing)
	for i := 0; i < g.N(); i++ {
		rand := ring.RandInt(big.NewInt(3)).Uint64()
		g.Set(i, rand)
	}
	return g
}

// NewRandomPolyVec constructs a vector, whose elements sampled from the given `sampler`.
func NewRandomPolyVec(size int, sampler PolySampler, baseRing *ring.Ring) *PolyVec {
	v := NewPolyVec(size, baseRing)
	v.Populate(func(i int) Poly {
		return *NewRandomPoly(sampler, baseRing)
	})
	return v
}

// NewRandomPolyMatrix constructs a 2D matrix, whose elements sampled from the given `sampler`.
func NewRandomPolyMatrix(rows int, cols int, sampler PolySampler, baseRing *ring.Ring) *PolyMatrix {
	A := NewPolyMatrix(rows, cols, baseRing)
	A.PopulateRows(func(_ int) PolyVec {
		return *NewRandomPolyVec(cols, sampler, baseRing)
	})
	return A
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

// NewRandomIntMatrix constructs a random 2D matrix of integers mod n.
func NewRandomIntMatrix(rows int, cols int, n *big.Int, baseRing *ring.Ring) *IntMatrix {
	A := NewIntMatrix(rows, cols, baseRing)
	A.PopulateRows(func(_ int) *IntVec {
		return NewRandomIntVec(cols, n, baseRing)
	})
	return A
}

// NewRandomIntMatrix constructs a random 2D matrix of integers mod n.
func NewRandomBinaryIntMatrix(rows int, cols int, baseRing *ring.Ring) *IntMatrix {
	A := NewIntMatrix(rows, cols, baseRing)
	A.PopulateRows(func(_ int) *IntVec {
		return NewRandomBinaryIntVec(cols, baseRing)
	})
	return A
}
