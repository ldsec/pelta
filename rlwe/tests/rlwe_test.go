package tests

import (
	"github.com/ldsec/codeBase/commitment/ens20"
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/ldsec/codeBase/commitment/math/rings"
	"github.com/ldsec/codeBase/commitment/rlwe"
	"testing"
)

func TestSaveLoadNTTTransform(t *testing.T) {
	settings := ens20.GetSimpleTestSettings()
	savedT, err := rlwe.SaveNTTTransform(settings.BaseRing, settings.Q, uint64(settings.LogD))
	if err != nil {
		t.Errorf(err.Error())
	}
	loadedT, err := rlwe.LoadNTTTransform(settings.BaseRing, settings.Q)
	if err != nil {
		t.Errorf(err.Error())
	}
	if !savedT.Eq(loadedT.MultiArray) {
		t.Errorf("saved and loaded transforms are not equal")
	}
}

func TestNTTTransform(t *testing.T) {
	settings := ens20.GetSimpleTestSettings()
	T, err := rlwe.LoadNTTTransform(settings.BaseRing, settings.Q)
	if err != nil {
		t.Errorf(err.Error())
	}
	// Create a test polynomial.
	poly := rings.NewZeroPolynomial(settings.BaseRing)
	for i := 0; i < poly.BaseRing.N; i++ {
		poly.SetCoefficient(i, uint64(i+1))
	}
	// Take the NTT transform by the transform matrix.
	polyNTT1 := rings.NewZIntVec(T.MulVec(poly.Coeffs(settings.Q).AsVec())).
		ToPoly(settings.BaseRing, settings.Q, true)
	polyNTT2 := poly.Copy().(rings.Polynomial).NTT()
	if !polyNTT1.Eq(polyNTT2) {
		t.Errorf("ntt transformation unsuccessful")
	}
}

func TestRLWETransformation(t *testing.T) {
	settings := ens20.GetSimpleTestSettings()
	s := math.NewRandomTernaryIntegerVector(settings.Kappa).
		ToPoly(settings.BaseRing, settings.Q, false)
	p1 := math.NewRandomPolynomial(settings.BaseRing, settings.UniformSampler)
	rlweParams := rlwe.NewRLWEParameters(settings)
	rlweProblem := rlwe.NewRLWEProblem(p1, s, settings.GaussianSampler, rlweParams)
	println("RLWE problem constructed")
	sisProblem := rlwe.RLWEToSIS(rlweProblem)
	println("A:", sisProblem.A.String())
}
