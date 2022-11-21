package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func verifySIS(t *testing.T, sis crypto.SISProblem) {
	u := sis.A.MulVec(sis.S)
	if !sis.U.Eq(u) {
		t.Errorf("SIS construction failed")
	}
}

func TestSISConstruction(t *testing.T) {
	config := crypto.GetDefaultConfig()
	n := 200
	m := 300
	A := fastmath.NewRandomIntMatrix(m, n, config.Q, config.BaseRing)
	s := fastmath.NewRandomIntVec(n, config.Q, config.BaseRing)
	sisProblem := crypto.NewSISProblem(A, s)
	verifySIS(t, sisProblem)
}

func TestSISSplit(t *testing.T) {
	config := crypto.GetDefaultConfig()
	n := 2 * config.D
	m := 2 * config.D
	A := fastmath.NewRandomIntMatrix(m, n, config.Q, config.BaseRing)
	s := fastmath.NewRandomIntVec(n, config.Q, config.BaseRing)
	sisProblem := crypto.NewSISProblem(A, s)
	smallRing := fastmath.ShortCommitmentRing(7)
	sisProblems := sisProblem.Split(config.RingParams, smallRing)
	for _, p := range sisProblems {
		verifySIS(t, p)
	}
}
