package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
)

func getBaseRing() *ring.Ring {
	return fastmath.BFVZeroLevelRing().BaseRing
}

func TestPolySumCoeffs(t *testing.T) {
	baseRing := getBaseRing()
	v := fastmath.NewPoly(baseRing)
	for i := 0; i < 5; i++ {
		v.SetForce(i, uint64(i+1))
	}
	actual := v.SumCoeffs()
	expected := uint64(1 + 2 + 3 + 4 + 5)
	if actual[0] != expected {
		t.Errorf("actual=%d, expected=%d", actual, expected)
	}
}

func TestIntVecRebaseLossless(t *testing.T) {
	largeRing := fastmath.BFVZeroLevelShortCommtRing(8)
	v := fastmath.NewIntVec(largeRing.D, largeRing.BaseRing)
	for i := 0; i < v.Size(); i++ {
		v.SetForce(i, uint64(i+1))
	}
	// Before rebase.
	// t.Logf(v.String())
	if len(v.UnderlyingPolys()) != 1 {
		t.Errorf("actual=%d, expected=%d", len(v.UnderlyingPolys()), 1)
	}
	if v.Size() != 256 {
		t.Errorf("actual=%d, expected=%d", v.Size(), 256)
	}
	v.RebaseLossless(fastmath.BFVZeroLevelShortCommtRing(4))
	// After rebase.
	if len(v.UnderlyingPolys()) != 16 {
		t.Errorf("actual=%d, expected=%d", len(v.UnderlyingPolys()), 16)
	}
	if v.Size() != 256 {
		t.Errorf("actual=%d, expected=%d", v.Size(), 256)
	}
	// t.Logf(v.String())
	for i := 0; i < v.Size(); i++ {
		if v.GetLevel(i, 0) != uint64(i+1) {
			t.Errorf("actual=%d, expected=%d", v.GetCoeff(i), i+1)
		}
	}
}

func TestIntVecDot(t *testing.T) {
	baseRing := getBaseRing()
	a := fastmath.NewIntVec(5, baseRing)
	b := fastmath.NewIntVec(5, baseRing)
	for i := 0; i < 5; i++ {
		a.SetForce(i, uint64(i+1))
		b.SetForce(i, uint64(i+2))
	}
	actual := a.Dot(b)
	expected := uint64(1*2 + 2*3 + 3*4 + 4*5 + 5*6)
	if actual[0] != expected {
		t.Errorf("actual=%d, expected=%d", actual, expected)
	}
}

func TestIntVecBigDot(t *testing.T) {
	baseRing := getBaseRing()
	size := baseRing.N*2 + 3
	a := fastmath.NewIntVec(size, baseRing)
	b := fastmath.NewIntVec(size, baseRing)
	expected := uint64(0)
	for i := 0; i < size; i++ {
		a.SetForce(i, uint64(i+1))
		b.SetForce(i, uint64(i+2))
		expected = (expected + uint64(i+1)*uint64(i+2)) % baseRing.ModulusAtLevel[0].Uint64()
	}
	actual := a.Dot(b)
	if actual[0] != expected {
		t.Errorf("actual=%d, expected=%d", actual, expected)
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
		t.Errorf("actualRows=%d, expectedRows=%d", at.Rows(), expectedRows)
	}
	if at.Cols() != expectedCols {
		t.Errorf("actualCols=%d, expectedCols=%d", at.Cols(), expectedCols)
	}
	expectedResult := [][]uint64{
		{0, 1, 2},
		{1, 2, 3},
		{2, 3, 4},
		{3, 4, 5},
		{4, 5, 6},
	}
	for i := 0; i < at.Rows(); i++ {
		for j := 0; j < at.Cols(); j++ {
			if at.GetLevel(i, j, 0) != expectedResult[i][j] {
				t.Errorf("actual[%d][%d] = %d, expected[%d][%d] = %d", i, j, at.GetCoeff(i, j), i, j, expectedResult[i][j])
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
	c := a.MulVec(v)
	expectedSize := 3
	if c.Size() != expectedSize {
		t.Errorf("actualSize=%d, expectedSize=%d", c.Size(), expectedSize)
	}
	expectedResult := []uint64{60, 80, 100}
	for i := 0; i < c.Size(); i++ {
		if c.GetLevel(i, 0) != expectedResult[i] {
			t.Errorf("actual[%d]=%d, expected[%d]=%d", i, c.GetCoeff(i), i, expectedResult[i])
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
	c := a.MulMat(b)
	expectedRows := 3
	expectedCols := 3
	if c.Rows() != expectedRows {
		t.Errorf("actualRows=%d, expectedRows=%d", c.Rows(), expectedRows)
	}
	if c.Cols() != expectedCols {
		t.Errorf("actualCols=%d, expectedCols=%d", c.Cols(), expectedCols)
	}
	expectedResult := [][]uint64{
		{0, 30, 60},
		{0, 40, 80},
		{0, 50, 100},
	}
	for i := 0; i < c.Rows(); i++ {
		for j := 0; j < c.Cols(); j++ {
			if c.GetLevel(i, j, 0) != expectedResult[i][j] {
				t.Errorf("actual[%d][%d]=%d, expected[%d][%d]=%d", i, j, c.GetCoeff(i, j), i, j, expectedResult[i][j])
			}
		}
	}
}

func TestIntMatrixRebaseRowsLossless(t *testing.T) {
	largeRing := fastmath.BFVZeroLevelShortCommtRing(8)
	m := fastmath.NewIntMatrix(largeRing.D, largeRing.D, largeRing.BaseRing)
	for i := 0; i < m.Rows(); i++ {
		for j := 0; j < m.Cols(); j++ {
			m.SetForce(i, j, uint64(i+1))
		}
	}
	// Before rebase.
	// t.Logf(m.String())
	for _, v := range m.RowsView() {
		if len(v.UnderlyingPolys()) != 1 {
			t.Errorf("actual=%d, expected=%d", len(v.UnderlyingPolys()), 1)
		}
		if v.Size() != 256 {
			t.Errorf("actual=%d, expected=%d", v.Size(), 256)
		}
	}
	m.RebaseRowsLossless(fastmath.BFVZeroLevelShortCommtRing(4))
	// After rebase.
	for _, v := range m.RowsView() {
		if len(v.UnderlyingPolys()) != 16 {
			t.Errorf("actual=%d, expected=%d", len(v.UnderlyingPolys()), 16)
		}
		if v.Size() != 256 {
			t.Errorf("actual=%d, expected=%d", v.Size(), 256)
		}
	}
	// t.Logf(m.String())
	for i := 0; i < m.Rows(); i++ {
		for j := 0; j < m.Cols(); j++ {
			if m.GetLevel(i, j, 0) != uint64(i+1) {
				t.Errorf("actual=%d, expected=%d", m.GetCoeff(i, j), i+1)
			}
		}
	}
}
