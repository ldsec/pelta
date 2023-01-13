package tests

import (
	"math"
	"math/big"
	"testing"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestComputeBasis(t *testing.T) {
	mod := uint64(5857)
	base := uint64(3)
	n := 8
	basis := fastmath.GenerateBasis(base, n, mod)
	for i := 0; i < n; i++ {
		expected := uint64(math.Pow(float64(base), float64(i)))
		if basis[i] != expected {
			t.Errorf("%d^%d != %d", base, i, expected)
		}
	}
}

func TestIntoBasis(t *testing.T) {
	base := uint64(3)
	num := uint64(174)
	logNum := 5
	basisRepr := fastmath.IntoBasisRepr(num, base, logNum)
	expected := []uint64{0, 1, 1, 0, 2}
	for i := 0; i < len(expected); i++ {
		if basisRepr[i] != expected[i] {
			t.Errorf("wrong basis returned!")
		}
	}
}

func TestTernaryDecomposition(t *testing.T) {
	baseRing := getBaseRing()
	vec := fastmath.NewIntVecFromSlice([]uint64{23, 174, 92, 3, 0}, baseRing)
	decomp, _ := fastmath.TernaryDecomposition(vec, big.NewInt(200), 5, baseRing)
	expectedDecomp := fastmath.NewIntMatrixFromSlice([][]uint64{
		{2, 1, 2, 0, 0},
		{0, 1, 1, 0, 2},
		{2, 0, 1, 0, 1},
		{0, 1, 0, 0, 0},
		{0, 0, 0, 0, 0},
	}, baseRing)
	if !decomp.Eq(expectedDecomp.Transposed()) {
		t.Errorf("wrong decomposition returned!")
	}
}
