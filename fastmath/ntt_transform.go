package fastmath

import (
	"fmt"
	"math"

	"github.com/ldsec/codeBase/commitment/logging"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// LoadNTTTransform gets the NTT transform associated with the given name. Either reads from the saved matrix file or generates anew.
func LoadNTTTransform(name string, baseRing *ring.Ring) *IntMatrix {
	return PersistentIntMatrix(name, func() *IntMatrix { return GenerateNTTTransform(baseRing) }, baseRing)
}

// GenerateNTTTransform computes and returns the integer NTT transformation matrix for the given base ring.
func GenerateNTTTransform(baseRing *ring.Ring) *IntMatrix {
	e := logging.LogExecStart("GenerateNTTTransform", "generation")
	defer e.LogExecEnd()
	Tlvls := make([][][]uint64, len(baseRing.ModulusAtLevel))
	for lvl := range Tlvls {
		Tlvls[lvl] = generateNTTTransformLevel(lvl, baseRing)
	}
	T := NewIntMatrix(baseRing.N, baseRing.N, baseRing)
	for i := 0; i < baseRing.N; i++ {
		for j := 0; j < baseRing.N; j++ {
			c := make([]uint64, len(baseRing.ModulusAtLevel))
			for lvl := range Tlvls {
				c[lvl] = Tlvls[lvl][i][j]
			}
			T.SetCoeff(i, j, c)
		}
	}
	return T
}

func generateNTTTransformLevel(lvl int, baseRing *ring.Ring) [][]uint64 {
	e := logging.LogExecStart("GenerateNTTTransformLevel", fmt.Sprintf("generation level %d", lvl))
	defer e.LogExecEnd()
	qUint := baseRing.Modulus[lvl]
	logN := int(math.Log2(float64(baseRing.N)))
	w := ring.InvMForm(baseRing.NttPsi[lvl][baseRing.N>>1], qUint, baseRing.MredParams[lvl])
	mask := uint64(2*baseRing.N - 1)
	T := make([][]uint64, baseRing.N)
	for i := 0; i < baseRing.N; i++ {
		T[i] = make([]uint64, baseRing.N)
		// Construct the transform row by row.
		twoirev := 2*utils.BitReverse64(uint64(i), uint64(logN)) + 1
		for j := 0; j < baseRing.N; j++ {
			gen := uint64(j) * twoirev & mask
			T[i][j] = ring.ModExp(w, gen, qUint)
		}
	}
	return T
}

// DiagMulMat performs the multiplication diag(v) * T in-place.
func (T *IntMatrix) DiagMulMat(v *IntVec) *IntMatrix {
	if T.Rows() != v.Size() {
		panic("incorrect sizes for diag(v) * T")
	}
	// Update the transformation matrix in-place without copying anything.
	for row := 0; row < T.Rows(); row++ {
		T.RowView(row).ScaleCoeff(v.GetCoeff(row))
	}
	return T
}
