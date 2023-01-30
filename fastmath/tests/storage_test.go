package tests

import (
	"os"
	"testing"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

func TestLoadSaveIntMatrix(t *testing.T) {
	testFileName := "test_matrix"
	baseRing := getBaseRing()
	matrix := fastmath.NewRandomIntMatrix(1000, 1000, baseRing.ModulusAtLevel[0], baseRing)
	// t.Logf("created matrix %s\n", matrix.String())
	err := fastmath.SaveIntMatrix(matrix, testFileName)
	if err != nil {
		t.Errorf("saving returned error %s", err.Error())
	}
	loadedMatrix, err := fastmath.LoadIntMatrix(testFileName, baseRing)
	if err != nil {
		t.Errorf("loading returned error %s", err.Error())
	}
	// t.Logf("loaded matrix: %s\n", loadedMatrix.String())
	// Make sure that the loaded matrix is correct.
	if !matrix.Eq(loadedMatrix) {
		t.Errorf("incorrect matrix")
	}
	// Clean up
	if os.Remove(testFileName) != nil {
		t.Logf("cleanup failed")
	}
}