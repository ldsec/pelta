package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/ring"
)

func getBaseRing() *ring.Ring {
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParams, _ := bfv.NewParametersFromLiteral(ringParamDef)
	return ringParams.CopyNew().Parameters.RingQP().RingQ
}

func TestPolySumCoeffs(t *testing.T) {
	baseRing := getBaseRing()
	v := fastmath.NewZeroPoly(baseRing)
	for i := 0; i < 5; i++ {
		v.Set(i, uint64(i+1))
	}
	actual := v.SumCoeffs(0)
	expected := uint64(1 + 2 + 3 + 4 + 5)
	if actual != expected {
		t.Errorf("TestPolySumCoeffs: actual=%d, expected=%d", actual, expected)
	}
}

func TestIntVecDot(t *testing.T) {
	baseRing := getBaseRing()
	a := fastmath.NewIntVec(5, baseRing)
	b := fastmath.NewIntVec(5, baseRing)
	for i := 0; i < 5; i++ {
		a.Set(i, uint64(i+1))
		b.Set(i, uint64(i+2))
	}
	actual := a.Dot(&b)
	expected := uint64(1*2 + 2*3 + 3*4 + 4*5 + 5*6)
	if actual != expected {
		t.Errorf("TestIntVecDot: actual=%d, expected=%d", actual, expected)
	}
}

func TestIntVecBigDot(t *testing.T) {
	baseRing := getBaseRing()
	size := baseRing.N*2 + 3
	a := fastmath.NewIntVec(size, baseRing)
	b := fastmath.NewIntVec(size, baseRing)
	expected := uint64(0)
	for i := 0; i < size; i++ {
		a.Set(i, uint64(i+1))
		b.Set(i, uint64(i+2))
		expected += uint64(i+1) * uint64(i+2)
	}
	actual := a.Dot(&b)
	if actual != expected {
		t.Errorf("TestIntVecBigDot: actual=%d, expected=%d", actual, expected)
	}
}

func TestIntMatrixTransposed(t *testing.T) {
	baseRing := getBaseRing()
	a := fastmath.NewIntMatrix(3, 5, baseRing)
	a.Populate(func(i, j int) uint64 {
		return uint64(i + j)
	})
	at := a.Transposed()
	expectedRows := 5
	expectedCols := 3
	if at.Rows() != expectedRows {
		t.Errorf("TestIntMatrixTranspose: actualRows=%d, expectedRows=%d", at.Rows(), expectedRows)
	}
	if at.Cols() != expectedCols {
		t.Errorf("TestIntMatrixTranspose: actualCols=%d, expectedCols=%d", at.Cols(), expectedCols)
	}
	expectedResult := [][]uint64{
		[]uint64{0, 1, 2},
		[]uint64{1, 2, 3},
		[]uint64{2, 3, 4},
		[]uint64{3, 4, 5},
		[]uint64{4, 5, 6},
	}
	for i := 0; i < at.Rows(); i++ {
		for j := 0; j < at.Cols(); j++ {
			if at.Get(i, j) != expectedResult[i][j] {
				t.Errorf("TestIntMatrixTranspose: actual[%d][%d] = %d, expected[%d][%d] = %d", i, j, at.Get(i, j), i, j, expectedResult[i][j])
			}
		}
	}
}

func TestIntMatrixMulVec(t *testing.T) {
	baseRing := getBaseRing()
	a := fastmath.NewIntMatrix(3, 5, baseRing)
	v := fastmath.NewIntVec(5, baseRing)
	a.Populate(func(i, j int) uint64 {
		return uint64(i + j)
	})
	v.Populate(func(i int) uint64 {
		return uint64(2 * i)
	})
	c := a.MulVec(&v)
	expectedSize := 3
	if c.Size() != expectedSize {
		t.Errorf("TestIntMatrixMulVec: actualSize=%d, expectedSize=%d", c.Size(), expectedSize)
	}
	expectedResult := []uint64{60, 80, 100}
	for i := 0; i < c.Size(); i++ {
		if c.Get(i) != expectedResult[i] {
			t.Errorf("TestIntMatrixMulVec: actual[%d]=%d, expected[%d]=%d", i, c.Get(i), i, expectedResult[i])
		}
	}
}

func TestIntMatrixMulMat(t *testing.T) {
	baseRing := getBaseRing()
	a := fastmath.NewIntMatrix(3, 5, baseRing)
	b := fastmath.NewIntMatrix(5, 3, baseRing)
	a.Populate(func(i, j int) uint64 {
		return uint64(i + j)
	})
	b.Populate(func(i, j int) uint64 {
		return uint64(i * j)
	})
	c := a.MulMat(&b)
	expectedRows := 3
	expectedCols := 3
	if c.Rows() != expectedRows {
		t.Errorf("TestIntMatrixMulMat: actualRows=%d, expectedRows=%d", c.Rows(), expectedRows)
	}
	if c.Cols() != expectedCols {
		t.Errorf("TestIntMatrixMulMat: actualCols=%d, expectedCols=%d", c.Cols(), expectedCols)
	}
	expectedResult := [][]uint64{
		[]uint64{0, 30, 60},
		[]uint64{0, 40, 80},
		[]uint64{0, 50, 100},
	}
	for i := 0; i < c.Rows(); i++ {
		for j := 0; j < c.Cols(); j++ {
			if c.Get(i, j) != expectedResult[i][j] {
				t.Errorf("TestIntMatrixMulMat: actual[%d][%d]=%d, expected[%d][%d]=%d", i, j, c.Get(i, j), i, j, expectedResult[i][j])
			}
		}
	}
}
