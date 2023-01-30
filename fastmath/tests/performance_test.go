package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestPerformanceRebase(t *testing.T) {
	bfvRing := fastmath.BFVFullRing()
	uni, _, _ := fastmath.GetSamplers(bfvRing, 128)
	A := fastmath.NewRandomIntMatrixFast(bfvRing.D, bfvRing.D, uni, bfvRing.BaseRing)

	smallRing := fastmath.BFVFullShortCommtRing(8)
	A.RebaseRowsLossless(smallRing)
}

func TestPerformanceMulVec(t *testing.T) {
	bfvRing := fastmath.BFVFullRing()
	uni, _, _ := fastmath.GetSamplers(bfvRing, 128)
	A := fastmath.NewRandomIntMatrixFast(bfvRing.D, bfvRing.D, uni, bfvRing.BaseRing)
	b := fastmath.NewRandomIntVecFast(bfvRing.D, uni, bfvRing.BaseRing)

	A.MulVec(b)
}
