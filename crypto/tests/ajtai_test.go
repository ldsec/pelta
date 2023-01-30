package tests

import (
	"math/big"
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func GetTestAjtaiConfig() crypto.AjtaiConfig {
	bfvRing := fastmath.BFVFullRing()
	p := big.NewInt(5857)
	return crypto.NewAjtaiConfig(p, bfvRing)
}

func TestAjtaiCommitments(t *testing.T) {
	config := GetTestAjtaiConfig()
	comSize := 4
	s := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	r := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	A := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	B := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	comQ := A.MulVec(s).Add(B.MulVec(r))
	comP := comQ.Copy().Reduce(config.P)
	comQActual, comPActual := crypto.GetAjtaiCommitments(A, B, s, r, config)
	if !comQActual.Eq(comQ) || !comPActual.Eq(comP) {
		t.Errorf("construction of com_p or com_q incorrect")
	}
}

func TestAjtaiKappa(t *testing.T) {
	config := GetTestAjtaiConfig()
	comSize := 4
	s := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	r := fastmath.NewRandomTernaryIntVec(config.D, config.BaseRing)
	A := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	B := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(A, B, s, r, config)
	kappa := crypto.GetAjtaiKappa(comP, comQ, config)
	// comP = comQ - kappa * p
	comPReconstructed := comQ.Copy().Add(kappa.Copy().
		Scale(config.P.Uint64()).
		Neg())
	if !comPReconstructed.Eq(comP) {
		t.Errorf("construction of kappa incorrect")
	}

}

func TestAjtaiEquationBuild(t *testing.T) {
	config := GetTestAjtaiConfig()
	uni, ter, _ := fastmath.GetSamplers(config.RingParams, 128)
	// Create another equation Mw = v.
	M := fastmath.NewRandomIntMatrixFast(config.D, config.D, uni, config.BaseRing)
	w := fastmath.NewRandomIntVec(config.D, config.Q, config.BaseRing)
	z := M.MulVec(w)
	// Create the Ajtai commitment.
	comSize := 4
	d := 256
	s := fastmath.NewRandomIntVecFast(d, ter, config.BaseRing)
	r := fastmath.NewRandomIntVecFast(d, ter, config.BaseRing)
	A := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	B := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(A, B, s, r, config)
	kappa := crypto.GetAjtaiKappa(comP, comQ, config)
	lrb := crypto.NewLinearRelationBuilder().
		AppendEqn(crypto.NewLinearEquation(z, w.Size()).AppendTerm(M, w)).
		AppendEqn(crypto.NewPaddedAjtaiEquation(comP, A, B, s, r, kappa, config))
	rel := lrb.Build(config.BaseRing)
	if !rel.IsValid() {
		t.Errorf("linearized ajtai eqn ill-formed")
	}
}

func TestAjtaiEquationBuildFast(t *testing.T) {
	config := GetTestAjtaiConfig()
	uni, ter, _ := fastmath.GetSamplers(config.RingParams, 128)
	// Create another equation Mw = v.
	M := fastmath.NewRandomIntMatrixFast(config.D, config.D, uni, config.BaseRing)
	w := fastmath.NewRandomIntVec(config.D, config.Q, config.BaseRing)
	z := M.MulVec(w)
	// Create the Ajtai commitment.
	comSize := 4
	d := 256
	s := fastmath.NewRandomIntVecFast(d, ter, config.BaseRing)
	r := fastmath.NewRandomIntVecFast(d, ter, config.BaseRing)
	A := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	B := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	comQ, comP := crypto.GetAjtaiCommitments(A, B, s, r, config)
	kappa := crypto.GetAjtaiKappa(comP, comQ, config)
	lrb := crypto.NewLinearRelationBuilder().
		AppendEqn(crypto.NewLinearEquation(z, w.Size()).AppendTerm(M, w)).
		AppendEqn(crypto.NewPaddedAjtaiEquation(comP, A, B, s, r, kappa, config))
	rel := lrb.BuildFast(config.BaseRing)
	if !rel.IsValid() {
		t.Errorf("linearized ajtai eqn ill-formed")
	}
}
