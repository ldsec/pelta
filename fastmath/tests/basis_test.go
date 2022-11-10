package tests

import (
	"math"
	"testing"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestComputeBasis(t *testing.T) {
	baseRing := getBaseRing()
	mod := baseRing.ModulusAtLevel[0].Uint64()
	base := uint64(3)
	n := 8
	basis := fastmath.GenerateBasis(base, n, mod, baseRing)
	for i := 0; i < n; i++ {
		expected := uint64(math.Pow(float64(base), float64(i)))
		if basis.Get(i) != expected {
			t.Errorf("TestComputeBasis: %d^%d != %d", base, i, expected)
		}
	}
}

func TestIntoBasis(t *testing.T) {
	baseRing := getBaseRing()
	base := uint64(3)
	num := uint64(174)
	logNum := 5
	basisRepr := fastmath.IntoBasisRepr(num, base, 200, logNum, baseRing)
	expected := fastmath.NewIntVecFromSlice([]uint64{0, 1, 1, 0, 2}, baseRing)
	if !basisRepr.Eq(&expected) {
		t.Errorf("TestIntoBasis: Wrong basis returned!")
	}
}

func TestTernaryDecomposition(t *testing.T) {
	baseRing := getBaseRing()
	mod := baseRing.ModulusAtLevel[0].Uint64()
	vec := fastmath.NewIntVecFromSlice([]uint64{23, 174, 92, 3, 0}, baseRing)
	decomp, _ := fastmath.TernaryDecomposition(vec, 200, 5, mod, baseRing)
	expectedDecomp := fastmath.NewIntMatrixFromSlice([][]uint64{
		{2, 1, 2, 0, 0},
		{0, 1, 1, 0, 2},
		{2, 0, 1, 0, 1},
		{0, 1, 0, 0, 0},
		{0, 0, 0, 0, 0},
	}, baseRing)
	if !decomp.Eq(&expectedDecomp) {
		t.Errorf("TestTernaryDecomposition: Wrong decomposition returned!")
	}
}
