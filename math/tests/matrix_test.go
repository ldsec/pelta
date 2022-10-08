package tests

import (
	"github.com/ldsec/codeBase/commitment/math"
	"math/big"
	"testing"
)

func TestMulVec(t *testing.T) {
	// ((1 2 3) (1 2 3) (1 2 3) (1 2 3) (1 2 3)) * (1 1 1) = (6 6 6 6 6)
	M := math.NewMatrixFromDimensions(5, 3).Populate(
		func(row int, col int) math.RingElement {
			return math.NewModInt((int64(col + 1)), big.NewInt(100))
		})
	x := math.NewVectorFromSize(3).Populate(
		func(i int) math.RingElement {
			return math.NewModInt((int64(1)), big.NewInt(100))
		})
	// Expected result
	y := math.NewVectorFromSize(5).Populate(
		func(i int) math.RingElement {
			return math.NewModInt((int64(6)), big.NewInt(100))
		})
	res := M.MulVec(x)
	for i := 0; i < len(y.Array); i++ {
		if y.Array[i].(*math.ModInt).Value.Cmp(&res.Array[i].(*math.ModInt).Value) != 0 {
			t.Errorf("Matrix.MulVec")
		}
	}
}

func TestTranspose(t *testing.T) {
	// ((1 1 1) (2 2 2) (3 3 3) (4 4 4))^T = ((1 2 3 4) (1 2 3 4) (1 2 3 4))
	M := math.NewMatrixFromDimensions(4, 3).Populate(
		func(row int, col int) math.RingElement {
			return math.NewModInt((int64(row + 1)), big.NewInt(100))
		})
	M.Transpose()
	if !(M.Rows() == 3 && M.Cols() == 4) {
		t.Errorf("Matrix.Transpose: wrong dimensions")
	}
	M.ForEach(func(el math.RingElement, row int, col int) {
		if el.(*math.ModInt).Value.Cmp(big.NewInt(int64(col)+1)) != 0 {
			t.Errorf("Matrix.Transpose: wrong values")
		}
	})
}
