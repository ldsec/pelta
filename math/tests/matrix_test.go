package tests

import (
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
	"math/big"
	"testing"
)

func TestMatrixDimensions(t *testing.T) {
	M := algebra.NewMatrixFromDimensions(5, 3).Populate(
		func(row int, col int) algebra.Element {
			return rings.NewZqInt(int64(col+1), big.NewInt(100))
		})
	if M.Rows() != 5 || M.Cols() != 3 {
		t.Errorf("Matrix.Dimensions")
	}
}

func TestMatrixMulVec(t *testing.T) {
	// ((1 2 3) (1 2 3) (1 2 3) (1 2 3) (1 2 3)) * (1 1 1) = (6 6 6 6 6)
	M := algebra.NewMatrixFromDimensions(5, 3).Populate(
		func(row int, col int) algebra.Element {
			return rings.NewZqInt(int64(col+1), big.NewInt(100))
		})
	x := algebra.NewVectorFromSize(3).Populate(
		func(i int) algebra.Element {
			return rings.NewZqInt(int64(1), big.NewInt(100))
		})
	// Expected result
	y := algebra.NewVectorFromSize(5).Populate(
		func(i int) algebra.Element {
			return rings.NewZqInt(int64(6), big.NewInt(100))
		})
	res := M.MulVec(x)
	for i := 0; i < len(y.Array); i++ {
		if y.Array[i].(*rings.ZInt).Value.Cmp(res.Array[i].(*rings.ZInt).Value) != 0 {
			t.Errorf("Matrix.MulVec")
		}
	}
}

func TestMatrixTranspose(t *testing.T) {
	// ((1 1 1) (2 2 2) (3 3 3) (4 4 4))^T = ((1 2 3 4) (1 2 3 4) (1 2 3 4))
	M := algebra.NewMatrixFromDimensions(4, 3).Populate(
		func(row int, col int) algebra.Element {
			return rings.NewZqInt(int64(row+1), big.NewInt(100))
		})
	M.Transpose()
	if !(M.Rows() == 3 && M.Cols() == 4) {
		t.Errorf("Matrix.Transpose: wrong dimensions")
	}
	M.ForEach(func(el algebra.Element, row int, col int) {
		if el.(*rings.ZInt).Value.Cmp(big.NewInt(int64(col)+1)) != 0 {
			t.Errorf("Matrix.Transpose: wrong values")
		}
	})
}
