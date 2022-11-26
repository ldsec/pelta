package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestAjtaiConstruction(t *testing.T) {
	config := crypto.GetDefaultConfig()
	comSize := 4
	s := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	r := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	aj := crypto.NewAjtaiCommitment(s, r, comSize, config)
	comQ := aj.A.MulVec(aj.S).
		Add(aj.B.MulVec(aj.R))
	if !aj.ComP.Eq(comQ.Copy().Reduce(aj.P)) {
		t.Errorf("construction of com_p incorrect")
	}
	// comP = comQ - kappa * p
	comP := comQ.Copy().Add(aj.Kappa.Copy().
		Scale(aj.P.Uint64()).
		Neg())
	if !aj.ComP.Eq(comP) {
		t.Errorf("construction of kappa incorrect")
	}
}

func TestAjtaiEmbedding(t *testing.T) {
	config := crypto.GetDefaultConfig()
	// Create the Ajtai commitment.
	comSize := 4
	s := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	r := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	aj := crypto.NewAjtaiCommitment(s, r, comSize, config)
	// Create the RLWE problem.
	p1 := fastmath.NewRandomPoly(config.UniformSampler, config.BaseRing)
	rlweParams := crypto.NewRLWEParameters(config.Q.Uint64(), config.D, uint64(config.Beta()), config.BaseRing)
	rlweProblem := crypto.NewRLWERelation(p1, s.UnderlyingPolys()[0].Copy(), config.GaussianSampler, rlweParams)
	// Convert into an SIS problem.
	sisProblem := crypto.RLWEToLinearRelation(rlweProblem)
	// Embed Ajtai into the problem.
	aj.EmbedIntoLinearRelation(&sisProblem, config.D, config.Q, config.BaseRing)
	// Make sure that the embedding results in a valid SIS problem.
	u := sisProblem.A.MulVec(sisProblem.S)
	if !sisProblem.U.Eq(u) {
		t.Errorf("embedded SIS has invalid construction")
	}
}
