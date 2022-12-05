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

func TestLRBSimple(t *testing.T) {

}