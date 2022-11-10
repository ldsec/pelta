package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestRLWEConstruction(t *testing.T) {
	config := crypto.GetDefaultConfig()
	// Create the RLWE problem.
	s := fastmath.NewRandomTernaryPoly(config.BaseRing)
	p1 := fastmath.NewRandomPoly(config.UniformSampler, config.BaseRing)
	rlweParams := crypto.NewRLWEParameters(config.Q.Uint64(), config.LogD, uint64(config.Beta()), config.BaseRing)
	rlweProblem := crypto.NewRLWEProblem(p1, s, config.GaussianSampler, rlweParams)
	// Check that p0 = -p1 * s + e
	p0 := rlweProblem.P1.Copy().NTT().Neg().Mul(rlweProblem.S.NTT()).Add(rlweProblem.E.NTT()).InvNTT()
	if !p0.Eq(rlweProblem.P0) {
		t.Errorf("RLWE construction is invalid")
	}
}

func TestRLWEErrorDecomposition(t *testing.T) {
	config := crypto.GetDefaultConfig()
	// Create the RLWE problem.
	s := fastmath.NewRandomTernaryPoly(config.BaseRing)
	p1 := fastmath.NewRandomPoly(config.UniformSampler, config.BaseRing)
	rlweParams := crypto.NewRLWEParameters(config.Q.Uint64(), config.LogD, uint64(config.Beta()), config.BaseRing)
	rlweProblem := crypto.NewRLWEProblem(p1, s, config.GaussianSampler, rlweParams)
	// Decompose error.
	e, b := rlweProblem.ErrorDecomposition()
	// Make sure that error = sum{e_i * b_i}
	reconstructedError := fastmath.NewIntVec(config.D, config.BaseRing)
	for i := 0; i < b.Size(); i++ {
		eRow := e.RowView(i).Copy()
		eRow.Scale(b.Get(i))
		reconstructedError.Add(eRow)
	}
	expectedError := rlweProblem.E.Coeffs()
	if !reconstructedError.Eq(expectedError) {
		t.Errorf("invalid decomposition")
	}
}

func TestRLWEToSIS(t *testing.T) {
	config := crypto.GetDefaultConfig()
	// Create the RLWE problem.
	s := fastmath.NewRandomTernaryPoly(config.BaseRing)
	p1 := fastmath.NewRandomPoly(config.UniformSampler, config.BaseRing)
	rlweParams := crypto.NewRLWEParameters(config.Q.Uint64(), config.LogD, uint64(config.Beta()), config.BaseRing)
	rlweProblem := crypto.NewRLWEProblem(p1, s, config.GaussianSampler, rlweParams)
	// Convert into an SIS problem.
	sisProblem := crypto.RLWEToSIS(rlweProblem)
	// Test that the equivalence p0 = -p1 * s + e holds.
	p0Coeffs1 := sisProblem.A.MulVec(sisProblem.S)
	p0 := rlweProblem.P0.Copy().NTT()
	p0Coeffs2 := p0.Coeffs()
	if !p0Coeffs1.Eq(p0Coeffs2) {
		t.Errorf("RLWE equivalence doesn't hold after the transformation")
	}
}
