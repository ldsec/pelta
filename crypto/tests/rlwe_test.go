package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func GetDefaultRLWEConfig() crypto.RLWEConfig {
	bfvRing := fastmath.BFVFullRingPN13()
	rlweErrorWidth := 128
	return crypto.NewRLWEConfig(rlweErrorWidth, bfvRing)
}

func TestRLWEP0(t *testing.T) {
	rlweConfig := GetDefaultRLWEConfig()
	uni, _, _ := fastmath.GetSamplers(rlweConfig.RingParams, uint64(rlweConfig.Delta))
	// Create the RLWE problem variables.
	p1 := fastmath.NewRandomPoly(uni, rlweConfig.BaseRing)
	s := fastmath.NewRandomTernaryPoly(rlweConfig.BaseRing)
	p0Actual, e := crypto.RLWESample(p1, s, rlweConfig)
	p0 := p1.Copy().NTT().Neg().Mul(s.Copy().NTT()).Add(e)
	if !p0Actual.Eq(p0) {
		t.Errorf("p_0 is incorrect")
	}
}

func TestRLWEErrorDecomposition(t *testing.T) {
	rlweParams := GetDefaultRLWEConfig()
	_, _, gau := fastmath.GetSamplers(rlweParams.RingParams, uint64(rlweParams.Delta))
	// Create a random error.
	err := fastmath.NewRandomPoly(gau, rlweParams.BaseRing)
	// Decompose error.
	e, b := crypto.RLWEErrorDecomposition(err, rlweParams)
	t.Logf("reconstructing")
	// Make sure that error = sum{e_i * b_i}
	reconstructedError := fastmath.NewIntVec(rlweParams.D, rlweParams.BaseRing)
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
	rlweParams := GetDefaultRLWEConfig()
	uni, _, _ := fastmath.GetSamplers(rlweParams.RingParams, uint64(rlweParams.Delta))
	// Create the RLWE problem variables.
	p1 := fastmath.NewRandomPoly(uni, rlweParams.BaseRing)
	s := fastmath.NewRandomTernaryPoly(rlweParams.BaseRing)
	p0, e := crypto.RLWESample(p1, s, rlweParams)
	// Construct the linear relation.
	T := fastmath.LoadNTTTransform("ntt.test", rlweParams.BaseRing)
	lrb := crypto.NewLinearRelationBuilder().AppendEqn(crypto.NewIndependentRLWE(p0, p1, s, e.InvNTT(), T, rlweParams))
	rlweRel := lrb.Build(rlweParams.BaseRing)
	verifyRelation(t, rlweRel)
}

func TestRLWEEquationBuildFast(t *testing.T) {
	rlweParams := GetDefaultRLWEConfig()
	uni, _, _ := fastmath.GetSamplers(rlweParams.RingParams, uint64(rlweParams.Delta))
	// Create the RLWE problem variables.
	p1 := fastmath.NewRandomPoly(uni, rlweParams.BaseRing)
	s := fastmath.NewRandomTernaryPoly(rlweParams.BaseRing)
	p0, e := crypto.RLWESample(p1, s, rlweParams)
	// Construct the linear relation.
	T := fastmath.LoadNTTTransform("ntt.test", rlweParams.BaseRing)
	lrb := crypto.NewLinearRelationBuilder().AppendEqn(crypto.NewIndependentRLWE(p0, p1, s, e.InvNTT(), T, rlweParams))
	rlweRel := lrb.BuildFast(rlweParams.BaseRing)
	verifyRelation(t, rlweRel)
}
