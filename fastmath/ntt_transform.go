package fastmath

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// GenerateNTTTransform computes and returns the integer NTT transformation matrix for the given base ring.
func GenerateNTTTransform(q uint64, logN int, baseRing *ring.Ring) IntMatrix {
	w := ring.InvMForm(baseRing.NttPsi[0][baseRing.N>>1], q, baseRing.MredParams[0])
	mask := uint64(2*baseRing.N - 1)
	T := NewIntMatrix(baseRing.N, baseRing.N, baseRing)
	T.PopulateRows(func(i int) IntVec {
		// Construct the transform row by row.
		twoirev := 2*utils.BitReverse64(uint64(i), uint64(logN)) + 1
		tRow := NewIntVec(baseRing.N, baseRing)
		tRow.Populate(func(j int) uint64 {
			gen := uint64(j) * twoirev & mask
			return ring.ModExp(w, gen, q)
		})
		return tRow
	})
	return T
}

func ExtendNTTTransform(t *IntMatrix, p *Poly) {
	// Update the transformation matrix in-place without copying anything.
	for row := 0; row < t.Rows(); row++ {
		tRow := t.RowView(row)
		tRow.Scale(p.Get(row, 0))
	}
}
