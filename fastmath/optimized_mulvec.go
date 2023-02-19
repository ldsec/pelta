package fastmath

import (
	"math/bits"
	"unsafe"

	"github.com/tuneinsight/lattigo/v4/ring"
)

func MulVecTransposeBlazingFast(A *IntMatrix, b *IntVec, baseRing *ring.Ring) *IntVec {
	out := NewIntVec(A.Cols(), baseRing)
	for lvl, qi := range baseRing.Modulus {
		ACoeffs := make([][]uint64, A.Rows())
		for i := 0; i < len(ACoeffs); i++ {
			ACoeffs[i] = A.RowView(i).GetWholeLevel(lvl)
		}
		bCoeffs := b.GetWholeLevel(lvl)
		outCoeffs := out.GetWholeLevel(lvl)
		MulVecTransposeBlazingFastLevel(ACoeffs, bCoeffs, outCoeffs, qi)
	}
	return out
}

func MulVecBlazingFast(A *IntMatrix, b *IntVec, baseRing *ring.Ring) *IntVec {
	out := NewIntVec(A.Rows(), baseRing)
	for lvl, qi := range baseRing.Modulus {
		ACoeffs := make([][]uint64, A.Rows())
		for i := 0; i < len(ACoeffs); i++ {
			ACoeffs[i] = A.RowView(i).GetWholeLevel(lvl)
		}
		bCoeffs := b.GetWholeLevel(lvl)
		outCoeffs := out.GetWholeLevel(lvl)
		for i := 0; i < A.Rows(); i++ {
			outCoeffs[i] = DotBlazingFastLevel(ACoeffs[i], bCoeffs, qi)
		}
	}
	return out
}

func DotBlazingFastLevel(V1, V2 []uint64, qi uint64) uint64 {
	mredConstant := ring.MRedParams(qi)
	bredConstant := ring.BRedParams(qi)
	out := uint64(0)
	for j := 0; j < len(V2); j += 1 {
		V2j := ring.MForm(V2[j], qi, bredConstant)
		out = ring.CRed(out+ring.MRed(V1[j], V2j, qi, mredConstant), qi)
	}
	return out
}

// MulVecTransposeBlazingFastLevel takes the matrix A and the vector V and
// returns U = A^T x V.
func MulVecTransposeBlazingFastLevel(A [][]uint64, V, U []uint64, qi uint64) {
	// U[j] = <A[*][j], V>
	mredConstant := ring.MRedParams(qi)
	bredConstant := ring.BRedParams(qi)
	rows := len(A)
	cols := len(A[0])
	//var wg sync.WaitGroup
	for i := 0; i < rows; i++ {
		Ai := A[i]
		Vi := ring.MForm(V[i], qi, bredConstant)
		//wg.Add(cols / 8)
		for j := 0; j < cols; j += 8 {
			func(j int) {
				x := (*[8]uint64)(unsafe.Pointer(&Ai[j]))
				z := (*[8]uint64)(unsafe.Pointer(&U[j]))
				// U[j] = U[j] + Ai[j] * v[i]
				z[0] = ring.CRed(z[0]+ring.MRed(x[0], Vi, qi, mredConstant), qi)
				z[1] = ring.CRed(z[1]+ring.MRed(x[1], Vi, qi, mredConstant), qi) // U[j+1]
				z[2] = ring.CRed(z[2]+ring.MRed(x[2], Vi, qi, mredConstant), qi) // U[j+2]
				z[3] = ring.CRed(z[3]+ring.MRed(x[3], Vi, qi, mredConstant), qi) // U[j+3]
				z[4] = ring.CRed(z[4]+ring.MRed(x[4], Vi, qi, mredConstant), qi) // ...
				z[5] = ring.CRed(z[5]+ring.MRed(x[5], Vi, qi, mredConstant), qi)
				z[6] = ring.CRed(z[6]+ring.MRed(x[6], Vi, qi, mredConstant), qi)
				z[7] = ring.CRed(z[7]+ring.MRed(x[7], Vi, qi, mredConstant), qi) // U[j+7]
				//wg.Done()
			}(j)
		}
		//wg.Wait()
	}
}

// MulATransposeByGammaV2 takes the matrix A and the vector Gamma and
// returns U = A^T x Gamma.
func MulATransposeByGammaV2(A [][]uint64, Gamma, U []uint64, qi uint64) {

	mredConstant := ring.MRedParams(qi) // MRedParams(qi) in previous Lattigo version
	bredConstant := ring.BRedParams(qi) // BRedParams(qi) in previous Lattigo version

	rows := len(A)
	cols := len(A[0])

	if rows != len(Gamma) || cols != len(U) {
		panic("dimensions A^T incompatible with Gamma or U")
	}

	clo := make([]uint64, cols)
	chi := make([]uint64, cols)

	log2qi := bits.Len64(qi)

	var maxiters = rows
	if 2*log2qi+bits.Len64(uint64(rows)) > 128 {
		maxiters = 1 << (128 - 2*log2qi)
	}

	mask := maxiters - 1

	for i := 0; i < rows; i++ {
		MulByScalarMontgomeryThenAddLazyU128(A[i], ring.MForm(Gamma[i], qi, bredConstant), clo, chi)

		if i&mask == maxiters-1 {
			ReduceU128ThenAdd(clo, chi, qi, mredConstant, bredConstant, U)
		}
	}

	if rows != maxiters-1 {
		ReduceU128ThenAdd(clo, chi, qi, mredConstant, bredConstant, U)
	}

}

// MulByScalarMontgomeryThenAddNoModU128 multiples the vector `a` with the scalar `scalar`
// and returns the lowers bits of the result in `clo` and the higher bits in `chi`.
// The scalar must be in the Mongtomery domain.
func MulByScalarMontgomeryThenAddLazyU128(a []uint64, scalar uint64, clo, chi []uint64) {

	N := len(a)

	var mhi, mlo, c uint64

	for i := 0; i < N; i += 8 {

		x := (*[8]uint64)(unsafe.Pointer(&a[i]))
		zlo := (*[8]uint64)(unsafe.Pointer(&clo[i]))
		zhi := (*[8]uint64)(unsafe.Pointer(&chi[i]))

		mhi, mlo = bits.Mul64(x[0], scalar)
		zlo[0], c = bits.Add64(zlo[0], mlo, 0)
		zhi[0] += mhi + c

		mhi, mlo = bits.Mul64(x[1], scalar)
		zlo[1], c = bits.Add64(zlo[1], mlo, 0)
		zhi[1] += mhi + c

		mhi, mlo = bits.Mul64(x[2], scalar)
		zlo[2], c = bits.Add64(zlo[2], mlo, 0)
		zhi[2] += mhi + c

		mhi, mlo = bits.Mul64(x[3], scalar)
		zlo[3], c = bits.Add64(zlo[3], mlo, 0)
		zhi[3] += mhi + c

		mhi, mlo = bits.Mul64(x[4], scalar)
		zlo[4], c = bits.Add64(zlo[4], mlo, 0)
		zhi[4] += mhi + c

		mhi, mlo = bits.Mul64(x[5], scalar)
		zlo[5], c = bits.Add64(zlo[5], mlo, 0)
		zhi[5] += mhi + c

		mhi, mlo = bits.Mul64(x[6], scalar)
		zlo[6], c = bits.Add64(zlo[6], mlo, 0)
		zhi[6] += mhi + c

		mhi, mlo = bits.Mul64(x[7], scalar)
		zlo[7], c = bits.Add64(zlo[7], mlo, 0)
		zhi[7] += mhi + c
	}
}

// ReduceU128ThenAdd takes the lower bits `alo` and higher bits `ahi` and
// finishes the modular reduction by `q`, adding the result on `b`.
func ReduceU128ThenAdd(alo, ahi []uint64, q, mredconstant uint64, bredconstant []uint64, b []uint64) {
	N := len(alo)

	var hhi uint64

	for i := 0; i < N; i += 8 {

		xlo := (*[8]uint64)(unsafe.Pointer(&alo[i]))
		xhi := (*[8]uint64)(unsafe.Pointer(&ahi[i]))
		z := (*[8]uint64)(unsafe.Pointer(&b[i]))

		hhi, _ = bits.Mul64(xlo[0]*mredconstant, q)
		z[0] = ring.BRedAdd(z[0]+xhi[0]-hhi+q, q, bredconstant)
		xlo[0], xhi[0] = 0, 0

		hhi, _ = bits.Mul64(xlo[1]*mredconstant, q)
		z[1] = ring.BRedAdd(z[1]+xhi[1]-hhi+q, q, bredconstant)
		xlo[1], xhi[1] = 0, 0

		hhi, _ = bits.Mul64(xlo[2]*mredconstant, q)
		z[2] = ring.BRedAdd(z[2]+xhi[2]-hhi+q, q, bredconstant)
		xlo[2], xhi[2] = 0, 0

		hhi, _ = bits.Mul64(xlo[3]*mredconstant, q)
		z[3] = ring.BRedAdd(z[3]+xhi[3]-hhi+q, q, bredconstant)
		xlo[3], xhi[3] = 0, 0

		hhi, _ = bits.Mul64(xlo[4]*mredconstant, q)
		z[4] = ring.BRedAdd(z[4]+xhi[4]-hhi+q, q, bredconstant)
		xlo[4], xhi[4] = 0, 0

		hhi, _ = bits.Mul64(xlo[5]*mredconstant, q)
		z[5] = ring.BRedAdd(z[5]+xhi[5]-hhi+q, q, bredconstant)
		xlo[5], xhi[5] = 0, 0

		hhi, _ = bits.Mul64(xlo[6]*mredconstant, q)
		z[6] = ring.BRedAdd(z[6]+xhi[6]-hhi+q, q, bredconstant)
		xlo[6], xhi[6] = 0, 0

		hhi, _ = bits.Mul64(xlo[7]*mredconstant, q)
		z[7] = ring.BRedAdd(z[7]+xhi[7]-hhi+q, q, bredconstant)
		xlo[7], xhi[7] = 0, 0
	}
}
