package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func verifyRelation(t *testing.T, rel *crypto.ImmutLinearRelation) {
	if !rel.IsValid() {
		t.Errorf("linear relation is ill-formed")
		// t.Logf(rel.U.String())
		// t.Logf(u.String())
	}
}

func TestLinRelConstruction(t *testing.T) {
	bfvRing := fastmath.BFVFullRingPN13()
	n := 200
	m := 300
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomIntVec(n, bfvRing.Q, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(fastmath.NewCachedIntMatrix(A), s).AsImmutable()
	verifyRelation(t, rel)
}

func TestLinRelRebase(t *testing.T) {
	largeRing := fastmath.BFVZeroLevelShortCommtRing(8)
	smallRing := fastmath.BFVZeroLevelShortCommtRing(4)
	n := largeRing.D
	m := largeRing.D
	A := fastmath.NewRandomIntMatrix(m, n, largeRing.Q, largeRing.BaseRing)
	s := fastmath.NewRandomIntVec(n, largeRing.Q, largeRing.BaseRing)
	rel := crypto.NewLinearRelation(fastmath.NewCachedIntMatrix(A), s).AsImmutable()
	t.Logf("testing over original ring...")
	verifyRelation(t, rel)
	rebasedRel := rel.Rebased(smallRing)
	t.Logf("testing over small ring...")
	verifyRelation(t, rebasedRel)
}

func TestLinRelAppendIndependent(t *testing.T) {
	bfvRing := fastmath.BFVFullRingPN13()
	n := 8
	m := 16
	// First relation: As = u
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomIntVec(n, bfvRing.Q, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(fastmath.NewCachedIntMatrix(A), s)
	verifyRelation(t, rel.AsImmutable())
	// Second relation: By = z
	np := 9
	mp := 13
	B := fastmath.NewRandomIntMatrix(mp, np, bfvRing.Q, bfvRing.BaseRing)
	y := fastmath.NewRandomIntVec(np, bfvRing.Q, bfvRing.BaseRing)
	rel2 := crypto.NewLinearRelation(fastmath.NewCachedIntMatrix(B), y)
	verifyRelation(t, rel2.AsImmutable())
	// Verify the appended version.
	rel.AppendIndependent(rel2)
	verifyRelation(t, rel.AsImmutable())
}

func TestLinRelExtend(t *testing.T) {
	bfvRing := fastmath.BFVFullRingPN13()
	n := 8
	m := 16
	// First relation: As = u
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomIntVec(n, bfvRing.Q, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(fastmath.NewCachedIntMatrix(A), s)
	verifyRelation(t, rel.AsImmutable())
	// Second relation: By = z
	np := 9
	mp := m
	B := fastmath.NewRandomIntMatrix(mp, np, bfvRing.Q, bfvRing.BaseRing)
	y := fastmath.NewRandomIntVec(np, bfvRing.Q, bfvRing.BaseRing)
	rel2 := crypto.NewLinearRelation(fastmath.NewCachedIntMatrix(B), y)
	verifyRelation(t, rel2.AsImmutable())
	// Verify the extended version.
	rel.Extend(rel2)
	verifyRelation(t, rel.AsImmutable())

}

func TestLinRelAppendDependentOnS(t *testing.T) {
	bfvRing := fastmath.BFVFullRingPN13()
	n := 8
	m := 16
	// First relation: As = u
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomIntVec(n, bfvRing.Q, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(fastmath.NewCachedIntMatrix(A), s)
	verifyRelation(t, rel.AsImmutable())
	// Second relation: Bs + y = z
	B := fastmath.NewRandomIntMatrix(n, n, bfvRing.Q, bfvRing.BaseRing)
	y := fastmath.NewRandomIntVec(n, bfvRing.Q, bfvRing.BaseRing)
	z := B.MulVec(s).Add(y)
	// Verify the appended version.
	rel.AppendDependentOnS(B, y, z)
	verifyRelation(t, rel.AsImmutable())
}
