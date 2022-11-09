package fastmath

import (
	"github.com/tuneinsight/lattigo/v4/ring"
)

// TernaryDecomposition computes an integer vector's ternary (0, 1, 2) decomposition.
// Returns the ternary decomposition matrix and the basis vector.
func TernaryDecomposition(v IntVec, logBeta int, mod uint64, baseRing *ring.Ring) (IntMatrix, IntVec) {
	base := uint64(3)
	basis := GenerateBasis(base, logBeta, mod, baseRing)
	decomposed := NewIntMatrix(v.Size(), logBeta, baseRing)
	decomposed.PopulateRows(func(i int) IntVec {
		basisRepr := IntoBasisRepr(v.Get(i), base, logBeta, baseRing)
		return basisRepr
	})
	return decomposed, basis
}

// IntoBasisRepr returns the given number's representation under the given base, where the digits are
// encoded in the returned vector in an LSB manner.
func IntoBasisRepr(num uint64, base uint64, numDigits int, baseRing *ring.Ring) IntVec {
	repr := make([]uint64, 0, numDigits)
	curr := num
	for i := 0; i < numDigits; i++ {
		rem := curr % base
		curr = curr / base
		repr = append(repr, rem)
	}
	return NewIntVecFromSlice(repr, baseRing)
}

// GenerateBasis returns a vector of multiplicants S.t. (base^0, base^1, ..., base^(n-1)) that can be used for
// base decomposition of vectors.
func GenerateBasis(base uint64, n int, mod uint64, baseRing *ring.Ring) IntVec {
	basis := NewIntVec(n, baseRing)
	basis.Populate(func(i int) uint64 {
		return ring.ModExp(base, uint64(i), mod)
	})
	return basis
}
