package tests

import (
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
	"github.com/ldsec/codeBase/commitment/rlwe"
	"math"
	"strconv"
	"testing"
)

func TestComputeBasis(t *testing.T) {
	base := rings.NewZInt(3)
	n := 8
	basis := rlwe.ComputeBasis(base, n)
	for i := 0; i < n; i++ {
		expected := int64(math.Pow(float64(base.Int64()), float64(i)))
		if basis.Element(i).(rings.ZInt).Int64() != expected {
			t.Errorf("Utils.ComputeBasis: %d^%d != %d", base, i, expected)
		}
	}
}

func TestIntoBasis(t *testing.T) {
	base := rings.NewZInt(3)
	num := rings.NewZInt(174)
	logNum := 5
	print(rlwe.IntoBasis(num, base, logNum).String())
	println(strconv.FormatInt(num.Int64(), int(base.Int64())))
}

func TestTernaryDecomposition(t *testing.T) {
	vec := rings.NewZIntVec(algebra.NewVectorFromSlice([]algebra.Element{
		rings.NewZInt(23),
		rings.NewZInt(174),
		rings.NewZInt(92),
		rings.NewZInt(3),
		rings.NewZInt(0),
	}))
	decomp, basis := rlwe.DecomposeIntoTernary(vec, 5)
	println("Decomp", decomp.String())
	println("Basis", basis.String())
}
