package optimizations

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

func OptimizedMulVec(m *fastmath.IntMatrix, v *fastmath.IntVec, mod *big.Int) *fastmath.IntVec {
	modUint := mod.Uint64()
	Ap := make([]uint64, m.Rows()*m.Cols())
	for i := 0; i < m.Rows(); i++ {
		for j := 0; j < m.Cols(); j++ {
			Ap[i*m.Cols()+j] = (m.Get(i, j) * v.Get(j)) % modUint
		}
	}
	Ab := make([]uint64, m.Rows())
	for i := 0; i < m.Rows(); i++ {
		for j := 0; j < m.Cols(); j++ {
			Ab[i] = (Ab[i] + Ap[i*m.Cols()+j]) % modUint
		}
	}
	out := fastmath.NewIntVec(m.Rows(), m.BaseRing())
	for i := 0; i < m.Rows(); i++ {
		out.Set(i, Ab[i])
	}
	return out
}

func OptimizedTranspose(m *fastmath.IntMatrix) *fastmath.IntMatrix {
	At := make([]uint64, m.Cols()*m.Rows())
	for row := 0; row < m.Rows(); row++ {
		for col := 0; col < m.Cols(); col++ {
			index := col*m.Rows() + row
			At[index] = m.Get(row, col)
		}
	}
	mt := fastmath.NewIntMatrix(m.Cols(), m.Rows(), m.BaseRing())
	for row := 0; row < mt.Rows(); row++ {
		for col := 0; col < mt.Cols(); col++ {
			mt.Set(row, col, At[row*mt.Cols()+col])
		}
	}
	return mt
} 