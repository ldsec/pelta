package tests

import (
	"math/big"
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestLinearEquationLinearize(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 2
	n := 3
	A := fastmath.NewRandomIntMatrix(m, n, big.NewInt(5), bfvRing.BaseRing)
	s := fastmath.NewRandomIntVec(n, big.NewInt(3), bfvRing.BaseRing)
	t.Logf(A.String())
	t.Logf(s.String())
	B := fastmath.NewRandomIntMatrix(m, n, big.NewInt(5), bfvRing.BaseRing)
	y := fastmath.NewRandomIntVec(n, big.NewInt(3), bfvRing.BaseRing)
	t.Logf(B.String())
	t.Logf(y.String())
	rel := crypto.NewLinearEquation(A.MulVec(s).Add(B.MulVec(y)), n).
		AppendTerm(A, s).
		AppendTerm(B, y).
		Linearize()
	t.Logf(rel.A.String())
	t.Logf(rel.S.String())
	t.Logf(rel.U.String())
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
	t.Logf(A.String())
	t.Logf(s.String())
	t.Logf(B.String())
	t.Logf(y.String())
	lrb.AppendEqn(crypto.NewLinearEquation(z, n).
		AppendTerm(A, s).
		AppendTerm(B, y))
	// Eqn 2: Cs + Dp = q
	C := fastmath.NewRandomIntMatrix(m, n, big.NewInt(5), bfvRing.BaseRing)
	D := fastmath.NewRandomIntMatrix(m, n, big.NewInt(5), bfvRing.BaseRing)
	p := fastmath.NewRandomIntVec(n, big.NewInt(3), bfvRing.BaseRing)
	q := C.MulVec(s).Add(D.MulVec(p))
	t.Logf(C.String())
	t.Logf(D.String())
	t.Logf(p.String())
	lrb.AppendEqn(crypto.NewLinearEquation(q, n).
		AppendDependentTerm(C, 0).
		AppendTerm(D, p))
	// Convert into linear relation (A || B || 0, C || 0 || D) (s y p) = (z q)
	linRel := lrb.Build(bfvRing.BaseRing)
	t.Logf(linRel.A.String())
	t.Logf(linRel.S.String())
	t.Logf(linRel.U.String())
	verifyRelation(t, linRel)
}