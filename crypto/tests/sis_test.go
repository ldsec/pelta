package tests

import (
	"math/big"
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestSISConstruction(t *testing.T) {
	config := crypto.GetDefaultConfig()
	n := big.NewInt(int64(config.D))
	A := fastmath.NewRandomIntMatrix(200, 200, n, config.BaseRing)
	s := fastmath.NewRandomIntVec(200, n, config.BaseRing)
	sisProblem := crypto.NewSISProblem(A, s)
	u := A.MulVec(s)
	if !sisProblem.U.Eq(u) {
		t.Errorf("SIS construction failed")
	}
}
