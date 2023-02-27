package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestPerformanceRebase(t *testing.T) {
	bfvRing := fastmath.BFVFullRingPN13()
	uni, _, _ := fastmath.GetSamplers(bfvRing, 128)
	A := fastmath.NewRandomIntMatrixFast(bfvRing.D, bfvRing.D, uni, bfvRing.BaseRing)

	smallRing := fastmath.BFVFullShortCommtRing(8)
	A.RebaseRowsLossless(smallRing)
}

func TestPerformanceMulVec(t *testing.T) {
	bfvRing := fastmath.BFVFullRingPN13()
	uni, _, _ := fastmath.GetSamplers(bfvRing, 128)
	numInstances := 2
	numRetries := 5
	A := make([]*fastmath.IntMatrix, numInstances)
	b := make([]*fastmath.IntVec, numInstances)
	for i := 0; i < numInstances; i++ {
		A[i] = fastmath.NewRandomIntMatrixFast(bfvRing.D, bfvRing.D, uni, bfvRing.BaseRing)
		b[i] = fastmath.NewRandomIntVecFast(bfvRing.D, uni, bfvRing.BaseRing)
	}
	for j := 0; j < numRetries; j++ {
		for i := 0; i < numInstances; i++ {
			A[i].MulVec(b[i])
		}
	}
}

// func TestPerformanceMulVecRebased(t *testing.T) {
// 	bfvRing := fastmath.BFVFullRingPN13()
// 	rebaseRing := fastmath.BFVFullShortCommtRing(8)
// 	uni, _, _ := fastmath.GetSamplers(bfvRing, 128)
// 	num := 2
// 	A := make([]*fastmath.IntMatrix, num)
// 	b := make([]*fastmath.IntVec, num)
// 	for i := 0; i < num; i++ {
// 		A[i] = fastmath.NewRandomIntMatrixFast(bfvRing.D, bfvRing.D, uni, bfvRing.BaseRing)
// 		b[i] = fastmath.NewRandomIntVecFast(bfvRing.D, uni, bfvRing.BaseRing)
// 	}
// 	for i := 0; i < num; i++ {
// 		A[i] = A[i].RebaseRowsLossless(rebaseRing).(*fastmath.IntMatrix)
// 		b[i] = b[i].RebaseLossless(rebaseRing)
// 	}
// 	for i := 0; i < num; i++ {
// 		A[i].MulVec(b[i])
// 	}
// }
