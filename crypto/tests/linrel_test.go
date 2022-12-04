package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func verifyRelation(t *testing.T, rel crypto.LinearRelation) {
	u := rel.A.MulVec(rel.S)
	if !rel.U.Eq(u) {
		t.Errorf("Linrel construction failed")
		t.Logf(rel.U.String())
		t.Logf(u.String())
	}
}

func TestLinRelConstruction(t *testing.T) {
	config := crypto.GetDefaultConfig()
	n := 200
	m := 300
	A := fastmath.NewRandomIntMatrix(m, n, config.Q, config.BaseRing)
	s := fastmath.NewRandomIntVec(n, config.Q, config.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	verifyRelation(t, rel)
}

func TestLinRelRebase(t *testing.T) {
	largeRing := fastmath.ShortCommitmentRing(8)
	smallRing := fastmath.ShortCommitmentRing(4)
	n := largeRing.D
	m := largeRing.D
	A := fastmath.NewRandomIntMatrix(m, n, largeRing.Q, largeRing.BaseRing)
	s := fastmath.NewRandomIntVec(n, largeRing.Q, largeRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	t.Logf("testing over original ring...")
	verifyRelation(t, rel)
	rebasedRel := rel.Rebase(smallRing)
	t.Logf("testing over small ring...")
	verifyRelation(t, rebasedRel)
}

func TestLinRelEmbedSecondary(t *testing.T) {
	config := crypto.GetDefaultConfig()
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
	rel.EmbedSecondary(B, y, z)
	verifyRelation(t, rel)
}