package fastmath

import (
	"github.com/tuneinsight/lattigo/v4/ring"
)

// TernaryDecomposition computes an integer vector's ternary (0, 1, 2) decomposition.
// Returns the ternary decomposition matrix and the basis vector.
func TernaryDecomposition(v *IntVec, logDelta int, baseRing *ring.Ring) (*IntMatrix, *IntVec) {
	base := uint64(3)
	basis := GenerateBasisCoeffs(base, logDelta, baseRing.Modulus)
	decomposed := NewIntMatrix(v.Size(), logDelta, baseRing)
	decomposed.PopulateRows(func(i int) *IntVec {
		// Decompose the ith coefficient of the error
		num := v.GetCoeff(i)
		basisRepr := IntoBasisReprCoeffs(num, base, logDelta)
		basisReprVec := NewIntVecFromCoeffSlice(basisRepr, baseRing)
		return basisReprVec
	})
	return decomposed.Transposed().(*IntMatrix), NewIntVecFromCoeffSlice(basis, baseRing)
}

// IntoBasisRepr returns the given number's representation under the given base, where the digits are
// encoded in the returned vector in an LSB manner.
func IntoBasisRepr(num uint64, base uint64, numDigits int) []uint64 {
	repr := make([]uint64, 0, numDigits)
	curr := num
	for i := 0; i < numDigits; i++ {
		rem := curr % base
		curr = curr / base
		repr = append(repr, rem)
	}
	return repr
}

// IntoBasisReprCoeffs performs the same operation with `IntoBasisRepr` but in parallel over the different levels
// of the given coefficient.
func IntoBasisReprCoeffs(num Coeff, base uint64, numDigits int) []Coeff {
	basisReprLvls := make([][]uint64, len(num))
	for lvl := 0; lvl < len(basisReprLvls); lvl++ {
		basisReprLvls[lvl] = IntoBasisRepr(num[lvl], base, numDigits)
	}
	coeffs := make([]Coeff, numDigits)
	for i := 0; i < len(coeffs); i++ {
		coeffs[i] = make([]uint64, len(num))
		for lvl := 0; lvl < len(basisReprLvls); lvl++ {
			coeffs[i][lvl] = basisReprLvls[lvl][i]
		}
	}
	return coeffs
}

// GenerateBasis returns a vector of multiplicants s.t. (base^0, base^1, ..., base^(n-1)) that can be used for
// base decomposition of vectors.
func GenerateBasis(base uint64, n int, mod uint64) []uint64 {
	basis := make([]uint64, n)
	for i := 0; i < n; i++ {
		basis[i] = ring.ModExp(base, uint64(i), mod)
	}
	return basis
}

// GenerateBasisCoeffs performs the same operation with `GenerateBasis` but in parallel over the different levels.
func GenerateBasisCoeffs(base uint64, n int, levelMods []uint64) []Coeff {
	basisLvl := make([][]uint64, len(levelMods))
	for lvl, q := range levelMods {
		basisLvl[lvl] = GenerateBasis(base, n, q)
	}
	coeffs := make([]Coeff, n)
	for i := 0; i < n; i++ {
		coeff := make([]uint64, len(levelMods))
		for lvl := 0; lvl < len(coeff); lvl++ {
			coeff[lvl] = basisLvl[lvl][i]
		}
		coeffs[i] = coeff
	}
	return coeffs
}
