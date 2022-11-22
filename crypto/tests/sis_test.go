package tests

import (
	"fmt"
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func verifySIS(t *testing.T, sis crypto.SISProblem) {
	u := sis.A.MulVec(sis.S)
	if !sis.U.Eq(u) {
		t.Errorf("SIS construction failed")
		fmt.Println(sis.U.String())
		fmt.Println(u.String())
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

func TestSISRebase(t *testing.T) {
	largeRing := fastmath.ShortCommitmentRing(8)
	smallRing := fastmath.ShortCommitmentRing(4)
	n := largeRing.D
	m := largeRing.D
	A := fastmath.NewRandomIntMatrix(m, n, largeRing.Q, largeRing.BaseRing)
	s := fastmath.NewRandomIntVec(n, largeRing.Q, largeRing.BaseRing)
	sisProblem := crypto.NewSISProblem(A, s)
	rebasedProblem := sisProblem.Rebase(smallRing)
	verifySIS(t, rebasedProblem)
}
