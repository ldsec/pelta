package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestAjtaiCommitments(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	params := crypto.NewAjtaiParameters(config.AjtaiMod, config.RingParams)
	comSize := 4
	s := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	r := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	A := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.AjtaiMod, config.BaseRing)
	B := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.AjtaiMod, config.BaseRing)
	comQ := A.MulVec(s).Add(B.MulVec(r))
	comP := comQ.Copy().Reduce(config.AjtaiMod)
	comQActual, comPActual := crypto.GetAjtaiCommitments(A, B, s, r, params)
	if !comQActual.Eq(comQ) || !comPActual.Eq(comP) {
		t.Errorf("construction of com_p or com_q incorrect")
	}
}

func TestAjtaiKappa(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	params := crypto.NewAjtaiParameters(config.AjtaiMod, config.RingParams)
	comSize := 4
	s := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	r := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	A := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.AjtaiMod, config.BaseRing)
	B := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.AjtaiMod, config.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(A, B, s, r, params)
	kappa := crypto.GetAjtaiKappa(comP, comQ, params)
	// comP = comQ - kappa * p
	comPReconstructed := comQ.Copy().Add(kappa.Copy().
		Scale(config.AjtaiMod.Uint64()).
		Neg())
	if !comPReconstructed.Eq(comP) {
		t.Errorf("construction of kappa incorrect")
	}

}

func TestAjtaiEquationBuild(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	params := crypto.NewAjtaiParameters(config.AjtaiMod, config.RingParams)
	// Create another equation Mw = v.
	M := fastmath.NewRandomIntMatrixFast(config.D, config.D, config.UniformSampler, config.BaseRing)
	w := fastmath.NewRandomIntVec(config.D, config.Q, config.BaseRing)
	z := M.MulVec(w)
	// Create the Ajtai commitment.
	comSize := 4
	d := 256
	s := fastmath.NewRandomIntVecFast(d, config.TernarySampler, config.BaseRing)
	r := fastmath.NewRandomIntVecFast(d, config.TernarySampler, config.BaseRing)
	A := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.AjtaiMod, config.BaseRing)
	B := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.AjtaiMod, config.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(A, B, s, r, params)
	kappa := crypto.GetAjtaiKappa(comP, comQ, params)
	lrb := crypto.NewLinearRelationBuilder().
		AppendEqn(crypto.NewLinearEquation(z, w.Size()).AppendTerm(M, w)).
		AppendEqn(crypto.NewPaddedAjtaiEquation(comP, A, B, s, r, kappa, params))
	rel := lrb.Build(config.BaseRing)
	if !rel.IsValid() {
		t.Errorf("linearized ajtai eqn ill-formed")
	}
}

func TestAjtaiEquationBuildFast(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	params := crypto.NewAjtaiParameters(config.AjtaiMod, config.RingParams)
	// Create another equation Mw = v.
	M := fastmath.NewRandomIntMatrixFast(config.D, config.D, config.UniformSampler, config.BaseRing)
	w := fastmath.NewRandomIntVec(config.D, config.Q, config.BaseRing)
	z := M.MulVec(w)
	// Create the Ajtai commitment.
	comSize := 4
	d := 256
	s := fastmath.NewRandomIntVecFast(d, config.TernarySampler, config.BaseRing)
	r := fastmath.NewRandomIntVecFast(d, config.TernarySampler, config.BaseRing)
	A := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.AjtaiMod, config.BaseRing)
	B := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.AjtaiMod, config.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(A, B, s, r, params)
	kappa := crypto.GetAjtaiKappa(comP, comQ, params)
	lrb := crypto.NewLinearRelationBuilder().
		AppendEqn(crypto.NewLinearEquation(z, w.Size()).AppendTerm(M, w)).
		AppendEqn(crypto.NewPaddedAjtaiEquation(comP, A, B, s, r, kappa, params))
	rel := lrb.BuildFast(config.BaseRing)
	if !rel.IsValid() {
		t.Errorf("linearized ajtai eqn ill-formed")
	}
}
