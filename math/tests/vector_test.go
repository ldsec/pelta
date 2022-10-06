package tests

import (
	"github.com/ldsec/codeBase/commitment/math"
	"math/big"
	"testing"
)

func TestDotProduct(t *testing.T) {
	// (1 2 3) * (1 2 3) = 14
	v1 := math.NewVectorFromSize(3).Populate(
		func(i int) math.RingElement {
			return math.NewModInt(big.NewInt(int64(i+1)), big.NewInt(100))
		})
	v2 := math.NewVectorFromSize(3).Populate(
		func(i int) math.RingElement {
			return math.NewModInt(big.NewInt(int64(i+1)), big.NewInt(100))
		})
	res := v1.DotProduct(v2).(*math.ModInt)
	if res.Value.Cmp(big.NewInt(14)) != 0 {
		t.Errorf("Vector.DotProduct")
	}
}
