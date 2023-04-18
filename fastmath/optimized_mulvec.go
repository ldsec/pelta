package fastmath

import (
	"github.com/ldsec/codeBase/commitment/misc"
	"github.com/ldsec/codeBase/commitment/relations"
	"math/bits"
	"unsafe"

	"github.com/tuneinsight/lattigo/v4/ring"
)

func MulVecTransposeBlazingFast(A *IntMatrix, b *IntVec, baseRing *ring.Ring) *IntVec {
	out := NewIntVec(A.Cols(), baseRing)
	// force unset the `unset` flag because we are directly working
	// with the underlying polys
	for _, p := range out.UnderlyingPolys() {
		p.unset = false
	}
	misc.ExecuteInParallel(len(baseRing.Modulus),
		func(lvl int) {
			qi := baseRing.Modulus[lvl]
			ACoeffs := make([][]uint64, A.Rows())
			for i := 0; i < len(ACoeffs); i++ {
				if A.RowView(i).IsUnset() {
					ACoeffs[i] = nil
				} else {
					ACoeffs[i] = A.RowView(i).GetWholeLevel(lvl)
				}
			}
			bCoeffs := b.GetWholeLevel(lvl)
			outCoeffs := make([]uint64, A.Cols())
			MulVecTransposeBlazingFastLevel(ACoeffs, bCoeffs, outCoeffs, qi)
			out.SetWholeLevel(lvl, outCoeffs)
		}, relations.Threads)
	return out
}

func MulVecBlazingFast(A *IntMatrix, b *IntVec, baseRing *ring.Ring) *IntVec {
	out := NewIntVec(A.Rows(), baseRing)
	for lvl, qi := range baseRing.Modulus {
		ACoeffs := make([][]uint64, A.Rows())
		bCoeffs := b.GetWholeLevel(lvl)
		for i := 0; i < len(ACoeffs); i++ {
			ACoeffs[i] = A.RowView(i).GetWholeLevel(lvl)
		}
		misc.ExecuteInParallel(A.Rows(), func(i int) {
			dotResult := DotBlazingFastLevel(ACoeffs[i], bCoeffs, qi)
			out.SetLevel(i, lvl, dotResult)
		}, relations.Threads)
	}
	return out
}

func DotBlazingFast(V1, V2 *IntVec, baseRing *ring.Ring) Coeff {
	out := NewZeroCoeff(len(baseRing.Modulus))
	misc.ExecuteInParallel(len(baseRing.Modulus), func(lvl int) {
		V1Coeffs := V1.GetWholeLevel(lvl)
		V2Coeffs := V2.GetWholeLevel(lvl)
		qi := baseRing.Modulus[lvl]
		out[lvl] = DotBlazingFastLevel(V1Coeffs, V2Coeffs, qi)
	}, relations.Threads)
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
	for i := 0; i < rows; i++ {
		Ai := A[i]
		if Ai == nil {
			continue
		}
		Vi := ring.MForm(V[i], qi, bredConstant)
		for j := 0; j < cols; j += 8 {
			Aij := (*[8]uint64)(unsafe.Pointer(&Ai[j]))
			Uj := (*[8]uint64)(unsafe.Pointer(&U[j]))
			// U[j] = U[j] + Ai[j] * v[i]
			Uj[0] = ring.CRed(Uj[0]+ring.MRed(Aij[0], Vi, qi, mredConstant), qi)
			Uj[1] = ring.CRed(Uj[1]+ring.MRed(Aij[1], Vi, qi, mredConstant), qi) // U[j+1]
			Uj[2] = ring.CRed(Uj[2]+ring.MRed(Aij[2], Vi, qi, mredConstant), qi) // U[j+2]
			Uj[3] = ring.CRed(Uj[3]+ring.MRed(Aij[3], Vi, qi, mredConstant), qi) // U[j+3]
			Uj[4] = ring.CRed(Uj[4]+ring.MRed(Aij[4], Vi, qi, mredConstant), qi) // ...
			Uj[5] = ring.CRed(Uj[5]+ring.MRed(Aij[5], Vi, qi, mredConstant), qi)
			Uj[6] = ring.CRed(Uj[6]+ring.MRed(Aij[6], Vi, qi, mredConstant), qi)
			Uj[7] = ring.CRed(Uj[7]+ring.MRed(Aij[7], Vi, qi, mredConstant), qi) // U[j+7]
		}
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
