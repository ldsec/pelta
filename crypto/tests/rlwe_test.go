package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestRLWEP0(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	// Create the RLWE problem variables.
	p1 := fastmath.NewRandomPoly(config.UniformSampler, config.BaseRing)
	s := fastmath.NewRandomTernaryPoly(config.BaseRing)
	e := fastmath.NewRandomPoly(config.GaussianSampler, config.BaseRing)
	p0 := p1.Copy().NTT().Neg().Mul(s.Copy().NTT()).Add(e.Copy().NTT())
	p0Actual := crypto.RLWESample(p1, s, e)
	if !p0Actual.Eq(p0) {
		t.Errorf("p_0 is incorrect")
	}
}

func TestRLWEErrorDecomposition(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	// Create a random error.
	err := fastmath.NewRandomPoly(config.GaussianSampler, config.BaseRing)
	rlweParams := crypto.NewRLWEParameters(config.RLWEErrorWidth, config.RingParams)
	// Decompose error.
	e, b := crypto.RLWEErrorDecomposition(err, rlweParams)
	t.Logf("reconstructing")
	// Make sure that error = sum{e_i * b_i}
	reconstructedError := fastmath.NewIntVec(config.D, config.BaseRing)
	for i := 0; i < b.Size(); i++ {
		eRow := e.RowView(i).Copy()
		eRow.ScaleCoeff(b.GetCoeff(i))
		reconstructedError.Add(eRow)
	}
	expectedError := err.Coeffs()
	if !reconstructedError.Eq(expectedError) {
		t.Errorf("invalid decomposition")
	}
}

func TestRLWEEquationBuild(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	// Create the RLWE problem variables.
	p1 := fastmath.NewRandomPoly(config.UniformSampler, config.BaseRing)
	s := fastmath.NewRandomTernaryPoly(config.BaseRing)
	e := fastmath.NewRandomPoly(config.GaussianSampler, config.BaseRing)
	p0 := crypto.RLWESample(p1, s, e)
	// Construct the linear relation.
	T := fastmath.LoadNTTTransform("ntt.test", config.BaseRing)
	rlweParams := crypto.NewRLWEParameters(config.RLWEErrorWidth, config.RingParams)
	lrb := crypto.NewLinearRelationBuilder().AppendEqn(crypto.NewIndependentRLWE(p0, p1, s, e, T, rlweParams))
	rlweRel := lrb.Build(config.BaseRing)
	verifyRelation(t, rlweRel)
}

func TestRLWEEquationBuildFast(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	// Create the RLWE problem variables.
	p1 := fastmath.NewRandomPoly(config.UniformSampler, config.BaseRing)
	s := fastmath.NewRandomTernaryPoly(config.BaseRing)
	e := fastmath.NewRandomPoly(config.GaussianSampler, config.BaseRing)
	p0 := crypto.RLWESample(p1, s, e)
	// Construct the linear relation.
	T := fastmath.LoadNTTTransform("ntt.test", config.BaseRing)
	rlweParams := crypto.NewRLWEParameters(config.RLWEErrorWidth, config.RingParams)
	lrb := crypto.NewLinearRelationBuilder().AppendEqn(crypto.NewIndependentRLWE(p0, p1, s, e, T, rlweParams))
	rlweRel := lrb.BuildFast(config.BaseRing)
	verifyRelation(t, rlweRel)
}
