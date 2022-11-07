package tests

import (
	"github.com/ldsec/codeBase/commitment/ens20"
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/ldsec/codeBase/commitment/math/rings"
	"github.com/ldsec/codeBase/commitment/rlwe"
	"testing"
)

func TestSaveLoadNTTTransform(t *testing.T) {
	// TODO fix
	settings := ens20.GetSimpleTestSettings()
	savedT, err := rlwe.SaveNTTTransform(settings.BaseRing, settings.Q, settings.LogD)
	if err != nil {
		t.Errorf(err.Error())
	}
	loadedT, err := rlwe.LoadNTTTransform(settings.BaseRing, settings.LogD)
	if err != nil {
		t.Errorf(err.Error())
	}
	if !savedT.Eq(loadedT.MultiArray) {
		t.Errorf("saved and loaded transforms are not equal")
	}
}

func TestNTTTransform(t *testing.T) {
	settings := ens20.GetSimpleTestSettings()
	T := rlwe.GenerateNTTTransform(settings.BaseRing, settings.Q, settings.LogD)
	// Create a test polynomial.
	poly := rings.NewZeroPolynomial(settings.BaseRing)
	for i := 0; i < poly.BaseRing.N; i++ {
		poly.SetCoefficient(i, uint64(i+1))
	}
	// Take the NTT transform by the transform matrix.
	polyNTT1 := T.Apply(poly)
	polyNTT2 := poly.Copy().(rings.Polynomial).NTT()
	if !polyNTT1.Eq(polyNTT2) {
		t.Errorf("ntt transformation unsuccessful")
	}
}

func TestNTTTransformScale(t *testing.T) {
	settings := ens20.GetSimpleTestSettings()
	T := rlwe.GenerateNTTTransform(settings.BaseRing, settings.Q, settings.LogD)
	// Create a test polynomial.
	poly := rings.NewZeroPolynomial(settings.BaseRing)
	for i := 0; i < poly.BaseRing.N; i++ {
		poly.SetCoefficient(i, uint64(i+1))
	}
	polyNTT1 := T.Scaled(5).Apply(poly)
	polyNTT2 := poly.Copy().Scale(5).(rings.Polynomial).NTT()
	if !polyNTT1.Eq(polyNTT2) {
		t.Errorf("ntt scaling unsuccessful")
	}
}

func TestNTTTransformExtension(t *testing.T) {
	settings := ens20.GetSimpleTestSettings()
	T := rlwe.GenerateNTTTransform(settings.BaseRing, settings.Q, settings.LogD)
	// Create a test polynomial.
	poly1 := rings.NewZeroPolynomial(settings.BaseRing)
	for i := 0; i < poly1.BaseRing.N; i++ {
		poly1.SetCoefficient(i, uint64(i+1))
	}
	poly2 := rings.NewZeroPolynomial(settings.BaseRing)
	for i := 0; i < poly2.BaseRing.N; i++ {
		poly2.SetCoefficient(i, uint64(i+1))
	}
	// Test the extensions.
	polyExtNTT1 := T.Extended(poly1).Apply(poly2).InvNTT()
	polyExtNTT2 := poly1.Copy().Mul(poly2)
	if !polyExtNTT1.Eq(polyExtNTT2) {
		t.Errorf("ntt extension unsuccessful")
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
