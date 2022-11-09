package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/ens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/rlwe"
)

func TestRLWEConstruction(t *testing.T) {
	settings := ens20.GetSimpleTestSettings()
	// Create the RLWE problem.
	s := fastmath.NewRandomTernaryPoly(settings.BaseRing)
	p1 := fastmath.NewRandomPoly(settings.UniformSampler, settings.BaseRing)
	rlweParams := rlwe.NewRLWEParameters(settings.Q.Uint64(), settings.LogD, uint64(settings.Beta()), settings.BaseRing)
	rlweProblem := rlwe.NewRLWEProblem(p1, s, settings.GaussianSampler, rlweParams)
	// Check that p0 = -p1 * s + e
	p0 := rlweProblem.P1.Copy()
	p0.NTT().Neg().MulCoeffs(rlweProblem.S.NTT()).Add(rlweProblem.E.NTT())
	p0.InvNTT()
	if !p0.Eq(&rlweProblem.P0) {
		t.Errorf("TestRLWEConstruction: RLWE construction is invalid!")
	}
}

func TestRLWEErrorDecomposition(t *testing.T) {
	settings := ens20.GetSimpleTestSettings()
	// Create the RLWE problem.
	s := fastmath.NewRandomTernaryPoly(settings.BaseRing)
	p1 := fastmath.NewRandomPoly(settings.UniformSampler, settings.BaseRing)
	rlweParams := rlwe.NewRLWEParameters(settings.Q.Uint64(), settings.LogD, uint64(settings.Beta()), settings.BaseRing)
	rlweProblem := rlwe.NewRLWEProblem(p1, s, settings.GaussianSampler, rlweParams)
	// Decompose error.
	e, b := rlweProblem.ErrorDecomposition()
	// Make sure that error = sum{e_i * b_i}
	reconstructedError := fastmath.NewIntVec(settings.D, settings.BaseRing)
	for i := 0; i < b.Size(); i++ {
		eRow := e.RowView(i).Copy()
		eRow.Scale(b.Get(i))
		reconstructedError.Add(&eRow)
	}
	expectedError := rlweProblem.E.Coeffs()
	if !reconstructedError.Eq(&expectedError) {
		t.Errorf("TestRLWEErrorDecomposition: Invalid decomposition")
	}
}

func TestRLWEToSIS(t *testing.T) {
	settings := ens20.GetSimpleTestSettings()
	// Create the RLWE problem.
	s := fastmath.NewRandomTernaryPoly(settings.BaseRing)
	p1 := fastmath.NewRandomPoly(settings.UniformSampler, settings.BaseRing)
	rlweParams := rlwe.NewRLWEParameters(settings.Q.Uint64(), settings.LogD, uint64(settings.Beta()), settings.BaseRing)
	rlweProblem := rlwe.NewRLWEProblem(p1, s, settings.GaussianSampler, rlweParams)
	// Convert into an SIS problem.
	sisProblem := rlwe.RLWEToSIS(rlweProblem)
	// Test that the equivalence p0 = -p1 * s + e holds.
	p0Coeffs1 := sisProblem.A.MulVec(&sisProblem.S)
	p0 := rlweProblem.P0.Copy()
	p0.NTT()
	p0Coeffs2 := p0.Coeffs()
	if !p0Coeffs1.Eq(&p0Coeffs2) {
		t.Errorf("TestRLWEToSIS: Equivalence doesn't hold after the transformation!")
	}
}
