package rlwe

import (
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
)

// DecomposeIntoTernary computes an integer vector v's ternary (0, 1, 2) decomposition.
// Returns the ternary decomposition matrix and the basis vector.
func DecomposeIntoTernary(v rings.ZIntVector, logBeta int) (algebra.Matrix, rings.ZIntVector) {
	base := rings.NewZInt(3)
	basis := ComputeBasis(base, logBeta)
	decomposed := algebra.NewMatrixFromDimensions(v.Length(), logBeta)
	v.ForEach(
		func(el algebra.Element, i int) {
			// Get the basis representation of the value.
			basisRepr := IntoBasis(el.(rings.ZInt), base, logBeta).AsVec()
			decomposed.SetRow(i, basisRepr)
		})
	return decomposed, basis
}

// IntoBasis returns the given number's representation under the given base, where the digits are
// encoded in the returned vector in an LSB manner.
func IntoBasis(num rings.ZInt, base rings.ZInt, logNum int) rings.ZIntVector {
	repr := make([]algebra.Element, 0, logNum)
	curr := num.Copy().(rings.ZInt)
	for i := 0; i < logNum; i++ {
		rem := curr.Copy().(rings.ZInt).Rem(base)
		curr.EuclideanDiv(base)
		repr = append(repr, rem)
	}
	return rings.NewZIntVec(algebra.NewVectorFromSlice(repr))
}

// ComputeBasis returns a vector of multiplicants s.t. (base^0, base^1, ..., base^(n-1)) that can be used for
// base decomposition of vectors.
func ComputeBasis(base rings.ZInt, n int) rings.ZIntVector {
	return rings.NewZIntVec(algebra.NewVectorFromSize(n).Populate(
		func(i int) algebra.Element {
			return base.Copy().Pow(uint64(i))
		}))
}
