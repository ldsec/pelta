package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestABPSimple(t *testing.T) {
	bfvRing := fastmath.BFVFullShortCommtRing(7)
	m := bfvRing.D
	n := bfvRing.D
	tau := bfvRing.D
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(n, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(rel, config)
	if !fastens20.ExecuteWithBoundProof(s, fastmath.NewSlice(0, n), tau, params) {
		t.Errorf("execution failed!")
	}
}

func TestABPSlice(t *testing.T) {
	bfvRing := fastmath.BFVFullShortCommtRing(7)
	m := bfvRing.D
	n := bfvRing.D * 2
	tau := bfvRing.D
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(n, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(rel, config)
	if !fastens20.ExecuteWithBoundProof(s, fastmath.NewSlice(0, n/2), tau, params) {
		t.Errorf("execution failed!")
	}
}
