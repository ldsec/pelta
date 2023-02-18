package tests

import (
	"math/big"
	"testing"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestPartitionedIntMatrixAsIntMatrix(t *testing.T) {
	testRing := fastmath.BFVFullShortCommtRing(7)
	uni, _, _ := fastmath.GetSamplers(testRing, 128)
	A1 := fastmath.NewRandomIntMatrixFast(testRing.D, testRing.D, uni, testRing.BaseRing)
	A2 := fastmath.NewRandomIntMatrixFast(testRing.D, testRing.D, uni, testRing.BaseRing)
	A3 := fastmath.NewRandomIntMatrixFast(testRing.D, testRing.D, uni, testRing.BaseRing)
	M := fastmath.NewEmptyPartitionedIntMatrix(testRing.BaseRing)
	M.Emplace(0, 0, A1)
	M.Emplace(1, 1, A2)
	M.Emplace(2, 2, A3)
	Mp := M.AsIntMatrix()
	if !Mp.SubsectionCopy(0, testRing.D, 0, testRing.D).Eq(A1) {
		t.Errorf("partitioned int matrix to concrete int matrix incorrect")
	}
	if !Mp.SubsectionCopy(testRing.D, 2*testRing.D, testRing.D, 2*testRing.D).Eq(A2) {
		t.Errorf("partitioned int matrix to concrete int matrix incorrect")
	}
	if !Mp.SubsectionCopy(2*testRing.D, 3*testRing.D, 2*testRing.D, 3*testRing.D).Eq(A3) {
		t.Errorf("partitioned int matrix to concrete int matrix incorrect")
	}
}

func TestPartitionedIntMatrixMulVec1(t *testing.T) {
	testRing := fastmath.BFVFullShortCommtRing(7)
	uni, _, _ := fastmath.GetSamplers(testRing, 128)
	A1 := fastmath.NewRandomIntMatrixFast(testRing.D, testRing.D, uni, testRing.BaseRing)
	A2 := fastmath.NewRandomIntMatrixFast(testRing.D, testRing.D, uni, testRing.BaseRing)
	A3 := fastmath.NewRandomIntMatrixFast(testRing.D, testRing.D, uni, testRing.BaseRing)
	M := fastmath.NewEmptyPartitionedIntMatrix(testRing.BaseRing)
	M.Emplace(0, 0, A1)
	M.Emplace(1, 1, A2)
	M.Emplace(2, 2, A3)
	v := fastmath.NewRandomIntVecFast(testRing.D*3, uni, testRing.BaseRing)
	r1 := M.MulVec(v)
	r2 := M.AsIntMatrix().MulVec(v)
	if !r1.Eq(r2) {
		t.Errorf("partitioned mulvec incorrect")
	}
}

func TestPartitionedIntMatrixMulVec2(t *testing.T) {
	testRing := fastmath.BFVFullShortCommtRing(7)
	uni, _, _ := fastmath.GetSamplers(testRing, 128)
	A1 := fastmath.NewRandomIntMatrixFast(testRing.D, testRing.D, uni, testRing.BaseRing)
	A2 := fastmath.NewRandomIntMatrixFast(testRing.D, testRing.D, uni, testRing.BaseRing)
	A3 := fastmath.NewRandomIntMatrixFast(testRing.D, testRing.D, uni, testRing.BaseRing)
	M := fastmath.NewEmptyPartitionedIntMatrix(testRing.BaseRing)
	M.Emplace(0, 0, A1)
	M.Emplace(1, 0, A2)
	M.Emplace(1, 1, A3)
	v := fastmath.NewRandomIntVecFast(testRing.D*2, uni, testRing.BaseRing)
	r1 := M.MulVec(v)
	r2 := M.AsIntMatrix().MulVec(v)
	if !r1.Eq(r2) {
		t.Errorf("partitioned mulvec incorrect")
	}
}

func TestPartitionedIntMatrixMulVec3(t *testing.T) {
	testRing := fastmath.BFVFullShortCommtRing(7)
	uni, _, _ := fastmath.GetSamplers(testRing, 128)
	A1 := fastmath.NewRandomIntMatrixFast(testRing.D, testRing.D, uni, testRing.BaseRing)
	A2 := fastmath.NewRandomIntMatrixFast(testRing.D, testRing.D, uni, testRing.BaseRing)
	A3 := fastmath.NewRandomIntMatrixFast(testRing.D, testRing.D, uni, testRing.BaseRing)
	M := fastmath.NewEmptyPartitionedIntMatrix(testRing.BaseRing)
	M.Emplace(0, 0, A1)
	M.Emplace(1, 0, A2)
	M.Emplace(1, 1, A3)
	A4 := fastmath.NewRandomIntMatrixFast(testRing.D, 2*testRing.D, uni, testRing.BaseRing)
	A5 := fastmath.NewRandomIntMatrixFast(testRing.D, testRing.D, uni, testRing.BaseRing)
	M2 := fastmath.NewEmptyPartitionedIntMatrix(testRing.BaseRing)
	M2.Emplace(0, 0, M)
	M2.Emplace(1, 0, A4)
	M2.Emplace(1, 1, A5)
	v := fastmath.NewRandomIntVecFast(M2.Cols(), uni, testRing.BaseRing)
	r1 := M2.MulVec(v)
	r2 := M2.AsIntMatrix().MulVec(v)
	if !r1.Eq(r2) {
		t.Errorf("partitioned mulvec incorrect")
	}
}

func TestPartitionedIntMatrixMulVecTranspose(t *testing.T) {
	testRing := fastmath.BFVFullShortCommtRing(7)
	//uni, _, _ := fastmath.GetSamplers(testRing, 128)
	A1 := fastmath.NewRandomIntMatrix(testRing.D, testRing.D, big.NewInt(5), testRing.BaseRing)
	A2 := fastmath.NewRandomIntMatrix(testRing.D, testRing.D, big.NewInt(5), testRing.BaseRing)
	M2 := fastmath.NewEmptyPartitionedIntMatrix(testRing.BaseRing)
	M2.Emplace(0, 0, A1)
	M2.Emplace(0, 1, A2)
	v := fastmath.NewRandomIntVec(M2.Rows(), big.NewInt(5), testRing.BaseRing)
	r1 := M2.MulVecTranspose(v)
	r2 := M2.AsIntMatrix().Transposed().MulVec(v)
	if !r1.Eq(r2) {
		t.Errorf("partitioned mulvec incorrect")
	}
}

func TestPartitionedIntMatrixMulVecTranspose2(t *testing.T) {
	testRing := fastmath.BFVFullShortCommtRing(7)
	//uni, _, _ := fastmath.GetSamplers(testRing, 128)
	A1 := fastmath.NewRandomIntMatrix(testRing.D, testRing.D, testRing.Q, testRing.BaseRing)
	A2 := fastmath.NewRandomIntMatrix(testRing.D, testRing.D, testRing.Q, testRing.BaseRing)
	A3 := fastmath.NewRandomIntMatrix(testRing.D, testRing.D, testRing.Q, testRing.BaseRing)
	A4 := fastmath.NewRandomIntMatrix(testRing.D, testRing.D, testRing.Q, testRing.BaseRing)
	A5 := fastmath.NewRandomIntMatrix(testRing.D, testRing.D, testRing.Q, testRing.BaseRing)
	M2 := fastmath.NewEmptyPartitionedIntMatrix(testRing.BaseRing)
	M2.Emplace(0, 0, A1)
	M2.Emplace(0, 1, A2)
	M2.Emplace(0, 2, A3)
	M2.Emplace(1, 3, A4)
	M2.Emplace(1, 4, A5)
	v := fastmath.NewRandomIntVec(M2.Rows(), big.NewInt(5), testRing.BaseRing)
	r1 := M2.MulVecTranspose(v)
	r2 := M2.AsIntMatrix().Transposed().MulVec(v)
	if !r1.Eq(r2) {
		t.Errorf("partitioned mulvec incorrect")
	}
}
