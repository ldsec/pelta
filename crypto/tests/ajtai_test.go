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
	A := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	B := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	aj := crypto.NewAjtaiCommitment(A, B, s, r, config.P, config.BaseRing)
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
	// config := crypto.GetDefaultConfig()
	// // Create the Ajtai commitment.
	// comSize := 4
	// s := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	// r := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	// A := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	// B := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	// aj := crypto.NewAjtaiCommitment(A, B, s, r, config.P, config.BaseRing)
	// // Create the RLWE problem.
	// p1 := fastmath.NewRandomPoly(config.UniformSampler, config.BaseRing)
	// err := fastmath.NewRandomPoly(config.GaussianSampler, config.BaseRing)
	// rlweParams := crypto.NewRLWEParameters(config.Q, config.D, uint64(config.Beta()), config.BaseRing)
	// rlweRel := crypto.NewRLWERelation(p1, s.UnderlyingPolys()[0].Copy(), err, rlweParams)
	// // Convert into an SIS problem.
	// T := fastmath.LoadNTTTransform("ntt_transform", config.Q, config.LogD, config.BaseRing)
	// e, b := rlweRel.ErrorDecomposition()
	// linRel := rlweRel.ToLinearRelation(e, b, T)
	// // Embed Ajtai into the problem.
	// aj.EmbedIntoLinearRelation(&linRel, config.D, config.Q, config.BaseRing)
	// // Make sure that the embedding results in a valid SIS problem.
	// u := linRel.A.MulVec(linRel.S)
	// if !linRel.U.Eq(u) {
	// 	t.Errorf("embedded SIS has invalid construction")
	// }
}
