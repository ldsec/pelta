package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestABPSimple(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 16
	n := 16
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(n, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	if !fastens20.ExecuteWithBoundProof(rel, 2) {
		t.Errorf("Protocol failed!")
	}
}
