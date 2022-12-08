package tests

import (
	"math/big"
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestLinEqnLinearize(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 2
	n := 3
	A := fastmath.NewRandomIntMatrix(m, n, big.NewInt(5), bfvRing.BaseRing)
	s := fastmath.NewRandomIntVec(n, big.NewInt(3), bfvRing.BaseRing)
	B := fastmath.NewRandomIntMatrix(m, n, big.NewInt(5), bfvRing.BaseRing)
	y := fastmath.NewRandomIntVec(n, big.NewInt(3), bfvRing.BaseRing)
	rel := crypto.NewLinearEquation(A.MulVec(s).Add(B.MulVec(y)), n).
		AppendTerm(A, s).
		AppendTerm(B, y).
		Linearize()
	verifyRelation(t, rel)
}

func TestLRBIndependentTerms(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 2
	n := 3
	lrb := crypto.NewLinearRelationBuilder()
	// Eqn 1: As + By = z
	A := fastmath.NewRandomIntMatrix(m, n, big.NewInt(5), bfvRing.BaseRing)
	s := fastmath.NewRandomIntVec(n, big.NewInt(3), bfvRing.BaseRing)
	B := fastmath.NewRandomIntMatrix(m, n, big.NewInt(5), bfvRing.BaseRing)
	y := fastmath.NewRandomIntVec(n, big.NewInt(3), bfvRing.BaseRing)
	z := A.MulVec(s).Add(B.MulVec(y))
	lrb.AppendEqn(crypto.NewLinearEquation(z, n).
		AppendTerm(A, s).
		AppendTerm(B, y))
	// Eqn 2: Cr + Dp = q
	C := fastmath.NewRandomIntMatrix(m, n, big.NewInt(5), bfvRing.BaseRing)
	r := fastmath.NewRandomIntVec(n, big.NewInt(3), bfvRing.BaseRing)
	D := fastmath.NewRandomIntMatrix(m, n, big.NewInt(5), bfvRing.BaseRing)
	p := fastmath.NewRandomIntVec(n, big.NewInt(3), bfvRing.BaseRing)
	q := C.MulVec(r).Add(D.MulVec(p))
	lrb.AppendEqn(crypto.NewLinearEquation(q, n).
		AppendTerm(C, r).
		AppendTerm(D, p))
	// Convert into linear relation
	linRel := lrb.Build(bfvRing.BaseRing)
	verifyRelation(t, linRel)
}

func TestLRBDependentTerms(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 2
	n := 2
	lrb := crypto.NewLinearRelationBuilder()
	// Eqn 1: As + By = z
	A := fastmath.NewRandomIntMatrix(m, n, big.NewInt(5), bfvRing.BaseRing)
	s := fastmath.NewRandomIntVec(n, big.NewInt(3), bfvRing.BaseRing)
	B := fastmath.NewRandomIntMatrix(m, n, big.NewInt(5), bfvRing.BaseRing)
	y := fastmath.NewRandomIntVec(n, big.NewInt(3), bfvRing.BaseRing)
	z := A.MulVec(s).Add(B.MulVec(y))
	lrb.AppendEqn(crypto.NewLinearEquation(z, n).
		AppendTerm(A, s).
		AppendTerm(B, y))
	// Eqn 2: Cs + Dp = q
	C := fastmath.NewRandomIntMatrix(m, n, big.NewInt(5), bfvRing.BaseRing)
	D := fastmath.NewRandomIntMatrix(m, n, big.NewInt(5), bfvRing.BaseRing)
	p := fastmath.NewRandomIntVec(n, big.NewInt(3), bfvRing.BaseRing)
	q := C.MulVec(s).Add(D.MulVec(p))
	lrb.AppendEqn(crypto.NewLinearEquation(q, n).
		AppendDependentTerm(C, 0).
		AppendTerm(D, p))
	// Convert into linear relation (A || B || 0, C || 0 || D) (s y p) = (z q)
	linRel := lrb.Build(bfvRing.BaseRing)
	verifyRelation(t, linRel)
}