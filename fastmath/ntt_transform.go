package fastmath

import (
	"errors"
	"fmt"
	"math/big"
	"os"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// LoadNTTTransform gets the NTT transform associated with the given name. Either reads from the saved matrix file or generates anew.
func LoadNTTTransform(name string, q *big.Int, logN int, baseRing *ring.Ring) *IntMatrix {
	var T *IntMatrix
	if _, err := os.Stat(name); errors.Is(err, os.ErrNotExist) {
		T = GenerateNTTTransform(q, logN, baseRing)
		err := SaveIntMatrix(T, name)
		if err != nil {
			fmt.Printf("couldn't save the ntt transform %s: %s\n", name, err.Error())
		}
	} else {
		T, err = LoadIntMatrix(name, baseRing)
		if err != nil {
			fmt.Printf("couldn't load the ntt transform %s, regenerating: %s\n", name, err.Error())
			T = GenerateNTTTransform(q, logN, baseRing)
		}
	}
	return T
}

// GenerateNTTTransform computes and returns the integer NTT transformation matrix for the given base ring.
func GenerateNTTTransform(q *big.Int, logN int, baseRing *ring.Ring) *IntMatrix {
	fmt.Printf("generating ntt transform... this may take a while")
	qUint := q.Uint64()
	w := ring.InvMForm(baseRing.NttPsi[0][baseRing.N>>1], qUint, baseRing.MredParams[0])
	mask := uint64(2*baseRing.N - 1)
	T := NewIntMatrix(baseRing.N, baseRing.N, baseRing)
	T.PopulateRows(func(i int) *IntVec {
		// Construct the transform row by row.
		twoirev := 2*utils.BitReverse64(uint64(i), uint64(logN)) + 1
		tRow := NewIntVec(baseRing.N, baseRing)
		tRow.Populate(func(j int) uint64 {
			gen := uint64(j) * twoirev & mask
			return ring.ModExp(w, gen, qUint)
		})
		return tRow
	})
	return T
}

// DiagMulMat performs the multiplication diag(v) * T.
func (T *IntMatrix) DiagMulMat(v *IntVec) *IntMatrix {
	if T.Rows() != v.Size() {
		panic("incorrect sizes for diag(v) * T")
	}
	// Update the transformation matrix in-place without copying anything.
	for row := 0; row < T.Rows(); row++ {
		T.RowView(row).Scale(v.Get(row))
	}
	return T
}
