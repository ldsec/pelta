package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestNTTTransformApply(t *testing.T) {
	baseRing := getBaseRing()
	T := fastmath.LoadNTTTransform("ntt.test", baseRing)
	// Create a test polynomial.
	poly := fastmath.NewPoly(baseRing)
	for i := 0; i < baseRing.N; i++ {
		poly.SetForce(i, uint64(i+1))
	}
	// Take the NTT transform by the transform matrix.
	polyCoeffs := poly.Coeffs()
	polyNTT1Coeffs := T.MulVec(polyCoeffs)
	// Take the NTT transform regularly.
	polyNTT2 := poly.Copy().NTT()
	polyNTT2Coeffs := polyNTT2.Coeffs()
	// Compare.
	if !polyNTT1Coeffs.Eq(polyNTT2Coeffs) {
		t.Errorf("NTT transformation unsuccessful")
		t.Errorf(polyNTT1Coeffs.String())
		t.Errorf(polyNTT2Coeffs.String())
	}
}

func TestNTTTransformScale(t *testing.T) {
	baseRing := getBaseRing()
	T := fastmath.LoadNTTTransform("ntt.test", baseRing)
	// T := fastmath.LoadNTTTransform("ntt.test", baseRing)
	// Create a test polynomial.
	poly := fastmath.NewPoly(baseRing)
	for i := 0; i < baseRing.N; i++ {
		poly.SetForce(i, uint64(i+1))
	}
	// Take the NTT transform by the scaled transform matrix.
	T.Scale(3)
	polyCoeffs := poly.Coeffs()
	polyNTT1Coeffs := T.MulVec(polyCoeffs)
	// Take the NTT transform regularly.
	polyNTT2 := poly.Copy().NTT()
	polyNTT2.Scale(3)
	polyNTT2Coeffs := polyNTT2.Coeffs()
	// Compare.
	if !polyNTT1Coeffs.Eq(polyNTT2Coeffs) {
		t.Errorf("NTT scaling unsuccessful")
	}
}

func TestNTTTransformExtend(t *testing.T) {
	baseRing := getBaseRing()
	T := fastmath.LoadNTTTransform("ntt.test", baseRing)
	// Create test polynomials.
	p0 := fastmath.NewPoly(baseRing)
	for i := 0; i < baseRing.N; i++ {
		p0.SetForce(i, uint64(i+1))
	}
	p1 := fastmath.NewPoly(baseRing)
	for i := 0; i < baseRing.N; i++ {
		p1.SetForce(i, uint64(i+1))
	}
	p0NTT := p0.Copy().NTT()
	// Take the NTT transform by the extended transform matrix.
	p1Coeffs := p1.Coeffs().Copy()
	p0p1NTTCoeffs1 := T.DiagMulMat(p0NTT.Coeffs()).MulVec(p1Coeffs)
	// Take the multiplication regularly.
	p0p1NTT := p0NTT.Copy().Mul(p1.NTT())
	p0p1NTTCoeffs2 := p0p1NTT.Coeffs()
	// Compare.
	if !p0p1NTTCoeffs1.Eq(p0p1NTTCoeffs2) {
		t.Logf(p0p1NTTCoeffs1.String())
		t.Logf(p0p1NTTCoeffs2.String())
		t.Errorf("NTT extension unsuccessful")
	}
}
