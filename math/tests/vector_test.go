package tests

import (
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
	"math/big"
	"testing"
)

func TestMatrixDotProduct(t *testing.T) {
	// (1 2 3) * (1 2 3) = 14
	v1 := algebra.NewVectorFromSize(3).Populate(
		func(i int) algebra.Element {
			return rings.NewModInt(int64(i+1), big.NewInt(100))
		})
	v2 := algebra.NewVectorFromSize(3).Populate(
		func(i int) algebra.Element {
			return rings.NewModInt(int64(i+1), big.NewInt(100))
		})
	res := v1.Dot(v2).(*rings.ZInt)
	if res.Value.Cmp(big.NewInt(14)) != 0 {
		t.Errorf("Vector.Dot")
	}
}
