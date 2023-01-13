package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestAjtaiCommitments(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	comSize := 4
	s := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	r := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	A := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	B := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	comQ := A.MulVec(s).Add(B.MulVec(r))
	comP := comQ.Copy().Reduce(config.P)
	comQActual, comPActual := crypto.GetAjtaiCommitments(A, B, s, r, config.P)
	if !comQActual.Eq(comQ) || !comPActual.Eq(comP) {
		t.Errorf("construction of com_p or com_q incorrect")
	}
}

func TestAjtaiKappa(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	comSize := 4
	s := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	r := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	A := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	B := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(A, B, s, r, config.P)
	kappa := crypto.GetAjtaiKappa(comP, comQ, config.P, config.BaseRing)
	// comP = comQ - kappa * p
	comPReconstructed := comQ.Copy().Add(kappa.Copy().
		Scale(config.P.Uint64()).
		Neg())
	if !comPReconstructed.Eq(comP) {
		t.Errorf("construction of kappa incorrect")
	}

}

func TestAjtaiEquation(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	// Create another equation Mw = v.
	M := fastmath.NewRandomIntMatrixFast(config.D, config.D, config.UniformSampler, config.BaseRing)
	w := fastmath.NewRandomIntVec(config.D, config.Q, config.BaseRing)
	z := M.MulVec(w)
	// Create the Ajtai commitment.
	comSize := 4
	d := 256
	s := fastmath.NewRandomIntVecFast(d, config.TernarySampler, config.BaseRing)
	r := fastmath.NewRandomIntVecFast(d, config.TernarySampler, config.BaseRing)
	A := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	B := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(A, B, s, r, config.P)
	kappa := crypto.GetAjtaiKappa(comP, comQ, config.P, config.BaseRing)
	lrb := crypto.NewLinearRelationBuilder().
		AppendEqn(crypto.NewLinearEquation(z, w.Size()).AppendTerm(M, w)).
		AppendEqn(crypto.NewPaddedAjtaiEquation(comP, A, B, s, r, kappa, config.P, config.BaseRing))
	rel := lrb.Build(config.BaseRing)
	if !rel.IsValid() {
		t.Errorf("linearized ajtai eqn ill-formed")
	}
}
