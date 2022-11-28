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
	e := fastmath.NewRandomPoly(config.GaussianSampler, config.BaseRing)
	rlweParams := crypto.NewRLWEParameters(config.Q, config.D, uint64(config.Beta()), config.BaseRing)
	rlweRel := crypto.NewRLWERelation(p1, s, e, rlweParams)
	// Check that p0 = -p1 * s + e
	p0 := rlweRel.P1.Copy().NTT().Neg().Mul(rlweRel.S.NTT()).Add(rlweRel.E.NTT()).InvNTT()
	if !p0.Eq(rlweRel.P0) {
		t.Errorf("RLWE construction is invalid")
	}
}

func TestRLWEErrorDecomposition(t *testing.T) {
	config := crypto.GetDefaultConfig()
	// Create the RLWE problem.
	s := fastmath.NewRandomTernaryPoly(config.BaseRing)
	p1 := fastmath.NewRandomPoly(config.UniformSampler, config.BaseRing)
	err := fastmath.NewRandomPoly(config.GaussianSampler, config.BaseRing)
	rlweParams := crypto.NewRLWEParameters(config.Q, config.D, uint64(config.Beta()), config.BaseRing)
	rlweRel := crypto.NewRLWERelation(p1, s, err, rlweParams)
	// Decompose error.
	e, b := rlweRel.ErrorDecomposition()
	// Make sure that error = sum{e_i * b_i}
	reconstructedError := fastmath.NewIntVec(config.D, config.BaseRing)
	for i := 0; i < b.Size(); i++ {
		eRow := e.RowView(i).Copy()
		eRow.Scale(b.Get(i))
		reconstructedError.Add(eRow)
	}
	expectedError := rlweRel.E.Coeffs()
	if !reconstructedError.Eq(expectedError) {
		t.Errorf("invalid decomposition")
	}
}

func TestRLWEToLinearRelation(t *testing.T) {
	config := crypto.GetDefaultConfig()
	// Create the RLWE problem.
	s := fastmath.NewRandomTernaryPoly(config.BaseRing)
	p1 := fastmath.NewRandomPoly(config.UniformSampler, config.BaseRing)
	err := fastmath.NewRandomPoly(config.GaussianSampler, config.BaseRing)
	rlweParams := crypto.NewRLWEParameters(config.Q, config.D, uint64(config.Beta()), config.BaseRing)
	rlweRel := crypto.NewRLWERelation(p1, s, err, rlweParams)
	// Convert into an SIS problem.
	T := fastmath.LoadNTTTransform("ntt_transform", config.Q, config.LogD, config.BaseRing)
	e, b := rlweRel.ErrorDecomposition()
	sisProblem := rlweRel.ToLinearRelation(e, b, T)
	// Test that the equivalence p0 = -p1 * s + e holds.
	p0Coeffs1 := sisProblem.A.MulVec(sisProblem.S)
	p0 := rlweRel.P0.Copy().NTT()
	p0Coeffs2 := p0.Coeffs()
	if !p0Coeffs1.Eq(p0Coeffs2) {
		t.Errorf("RLWE equivalence doesn't hold after the transformation")
	}
}
