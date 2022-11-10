package tests

import (
	"math"
	"testing"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestNTTTransform(t *testing.T) {
	baseRing := getBaseRing()
	q := baseRing.ModulusAtLevel[0].Uint64()
	logD := int(math.Log2(float64(baseRing.N)))
	T := fastmath.GenerateNTTTransform(q, logD, baseRing)
	// Create a test polynomial.
	poly := fastmath.NewZeroPoly(baseRing)
	for i := 0; i < baseRing.N; i++ {
		poly.Set(i, uint64(i+1))
	}
	// Take the NTT transform by the transform matrix.
	polyCoeffs := poly.Coeffs()
	polyNTT1Coeffs := T.MulVec(&polyCoeffs)
	// Take the NTT transform regularly.
	polyNTT2 := poly.Copy()
	polyNTT2.NTT()
	polyNTT2Coeffs := polyNTT2.Coeffs()
	// Compare.
	if !polyNTT1Coeffs.Eq(&polyNTT2Coeffs) {
		t.Errorf("TestNTTTransform: NTT transformation unsuccessful")
	}
}

func TestNTTTransformScale(t *testing.T) {
	baseRing := getBaseRing()
	q := baseRing.ModulusAtLevel[0].Uint64()
	logD := int(math.Log2(float64(baseRing.N)))
	T := fastmath.GenerateNTTTransform(q, logD, baseRing)
	// Create a test polynomial.
	poly := fastmath.NewZeroPoly(baseRing)
	for i := 0; i < baseRing.N; i++ {
		poly.Set(i, uint64(i+1))
	}
	// Take the NTT transform by the scaled transform matrix.
	T.Scale(3)
	polyCoeffs := poly.Coeffs()
	polyNTT1Coeffs := T.MulVec(&polyCoeffs)
	// Take the NTT transform regularly.
	polyNTT2 := poly.Copy()
	polyNTT2.Scale(3)
	polyNTT2.NTT()
	polyNTT2Coeffs := polyNTT2.Coeffs()
	// Compare.
	if !polyNTT1Coeffs.Eq(&polyNTT2Coeffs) {
		t.Errorf("TestNTTTransformScale: NTT scaling unsuccessful")
	}
}

func TestNTTTransformExtend(t *testing.T) {
	baseRing := getBaseRing()
	q := baseRing.ModulusAtLevel[0].Uint64()
	logD := int(math.Log2(float64(baseRing.N)))
	T := fastmath.GenerateNTTTransform(q, logD, baseRing)
	// Create test polynomials.
	p0 := fastmath.NewZeroPoly(baseRing)
	for i := 0; i < baseRing.N; i++ {
		p0.Set(i, uint64(i+1))
	}
	p1 := fastmath.NewZeroPoly(baseRing)
	for i := 0; i < baseRing.N; i++ {
		p1.Set(i, uint64(i+1))
	}
	p0.NTT()
	// Take the NTT transform by the extended transform matrix.
	fastmath.ExtendNTTTransform(&T, &p0)
	p1Coeffs := p1.Coeffs()
	p0p1NTTCoeffs1 := T.MulVec(&p1Coeffs)
	// Take the NTT transform regularly.
	p0p1NTT := p0.Copy().MulCoeffs(&p1)
	p0p1NTTCoeffs2 := p0p1NTT.Coeffs()
	// Compare.
	if !p0p1NTTCoeffs1.Eq(&p0p1NTTCoeffs2) {
		t.Errorf("TestNTTTransformExtend: NTT extension unsuccessful")
	}
}
