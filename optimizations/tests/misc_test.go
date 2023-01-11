package tests

import (
	"errors"
	"fmt"
	"os"
	"testing"
	"time"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/optimizations"
)

var matrix_name string = "random_matrix"
var vec_name string = "random_vec"
var m int = 7 * 8192
var n int = 2 * 8192

func generateTestObjects(config crypto.CryptoConfig) error {
	fmt.Println("generating test objects...")
	uni, _, _ := crypto.GetSamplers(config.RingParams, config.Delta1)
	A := fastmath.NewRandomIntMatrixFast(m, n, uni, config.BaseRing)
	b := fastmath.NewRandomIntVecFast(m, uni, config.BaseRing)
	err := fastmath.SaveIntMatrix(A, matrix_name)
	if err != nil {
		return err
	}
	err = fastmath.SaveIntVec(b, vec_name)
	if err != nil {
		return err
	}
	fmt.Println("done")
	return nil
}

func getRandomMatrix() *fastmath.IntMatrix {
	config := crypto.GetDefaultCryptoConfig()
	if _, err := os.Stat(matrix_name); errors.Is(err, os.ErrNotExist) {
		generateTestObjects(config)
	}
	fmt.Println("loading matrix")
	m, _ := fastmath.LoadIntMatrix(matrix_name, config.BaseRing)
	return m
}

func getRandomVec() *fastmath.IntVec {
	config := crypto.GetDefaultCryptoConfig()
	if _, err := os.Stat(vec_name); errors.Is(err, os.ErrNotExist) {
		generateTestObjects(config)
	}
	fmt.Println("loading vector")
	m, _ := fastmath.LoadIntVec(vec_name, config.BaseRing)
	return m
}

func TestOptimizedTranspose(t *testing.T) {
	A := getRandomMatrix()
	var AtOpt, AtUnopt *fastmath.IntMatrix
	{
		t.Logf("* unoptimized")
		t0 := time.Now()
		AtUnopt = A.Transposed()
		t1 := time.Now()
		t.Logf("took %dms", t1.Sub(t0).Milliseconds())
	}
	{
		t.Logf("* optimized")
		t0 := time.Now()
		AtOpt = optimizations.OptimizedTranspose(A)
		t1 := time.Now()
		t.Logf("took %dms", t1.Sub(t0).Milliseconds())
	}
	if !AtOpt.Eq(AtUnopt) {
		t.Logf("not equal")
	}
}

func TestOptimizedMulVec(t *testing.T) {
	A := optimizations.OptimizedTranspose(getRandomMatrix())
	b := getRandomVec()
	var AbOpt, AbUnopt *fastmath.IntVec
	{
		t.Logf("* unoptimized")
		t0 := time.Now()
		AbUnopt = A.MulVec(b)
		t1 := time.Now()
		t.Logf("took %dms", t1.Sub(t0).Milliseconds())
	}
	{
		t.Logf("* optimized")
		t0 := time.Now()
		AbOpt = optimizations.OptimizedMulVec(A, b, A.BaseRing().ModulusAtLevel[0])
		t1 := time.Now()
		t.Logf("took %dms", t1.Sub(t0).Milliseconds())
	}
	if !AbOpt.Eq(AbUnopt) {
		t.Logf("not equal")
	}
}
