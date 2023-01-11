package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func verifyRelation(t *testing.T, rel crypto.LinearRelation) {
	if !rel.IsValid() {
		t.Errorf("linear relation is ill-formed")
		// t.Logf(rel.U.String())
		// t.Logf(u.String())
	}
}

func TestLinRelConstruction(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	n := 200
	m := 300
	A := fastmath.NewRandomIntMatrix(m, n, config.Q, config.BaseRing)
	s := fastmath.NewRandomIntVec(n, config.Q, config.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	verifyRelation(t, rel)
}

func TestLinRelRebase(t *testing.T) {
	largeRing := fastmath.BFVZeroLevelShortCommtRing(8)
	smallRing := fastmath.BFVZeroLevelShortCommtRing(4)
	n := largeRing.D
	m := largeRing.D
	A := fastmath.NewRandomIntMatrix(m, n, largeRing.Q, largeRing.BaseRing)
	s := fastmath.NewRandomIntVec(n, largeRing.Q, largeRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	t.Logf("testing over original ring...")
	verifyRelation(t, rel)
	rebasedRel := rel.Rebased(smallRing)
	t.Logf("testing over small ring...")
	verifyRelation(t, rebasedRel)
}

func TestLinRelAppendIndependent(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	n := 8
	m := 16
	// First relation: As = u
	A := fastmath.NewRandomIntMatrix(m, n, config.Q, config.BaseRing)
	s := fastmath.NewRandomIntVec(n, config.Q, config.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	verifyRelation(t, rel)
	// Second relation: By = z
	np := 9
	mp := 13
	B := fastmath.NewRandomIntMatrix(mp, np, config.Q, config.BaseRing)
	y := fastmath.NewRandomIntVec(np, config.Q, config.BaseRing)
	rel2 := crypto.NewLinearRelation(B, y)
	verifyRelation(t, rel2)
	// Verify the appended version.
	rel.AppendIndependent(rel2)
	verifyRelation(t, rel)
}

func TestLinRelExtend(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	n := 8
	m := 16
	// First relation: As = u
	A := fastmath.NewRandomIntMatrix(m, n, config.Q, config.BaseRing)
	s := fastmath.NewRandomIntVec(n, config.Q, config.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	verifyRelation(t, rel)
	// Second relation: By = z
	np := 9
	mp := m
	B := fastmath.NewRandomIntMatrix(mp, np, config.Q, config.BaseRing)
	y := fastmath.NewRandomIntVec(np, config.Q, config.BaseRing)
	rel2 := crypto.NewLinearRelation(B, y)
	verifyRelation(t, rel2)
	// Verify the extended version.
	rel.Extend(rel2)
	verifyRelation(t, rel)

}

func TestLinRelAppendDependentOnS(t *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	n := 8
	m := 16
	// First relation: As = u
	A := fastmath.NewRandomIntMatrix(m, n, config.Q, config.BaseRing)
	s := fastmath.NewRandomIntVec(n, config.Q, config.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	verifyRelation(t, rel)
	// Second relation: Bs + y = z
	B := fastmath.NewRandomIntMatrix(n, n, config.Q, config.BaseRing)
	y := fastmath.NewRandomIntVec(n, config.Q, config.BaseRing)
	z := B.MulVec(s).Add(y)
	// Verify the appended version.
	rel.AppendDependentOnS(B, y, z)
	verifyRelation(t, rel)
}
