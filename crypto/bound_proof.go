package crypto

import (
	"fmt"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
)

func CreateABPChallenge(tau, l int, ternarySampler fastmath.PolySampler, baseRing *ring.Ring) *fastmath.PolyNTTMatrix {
	d := baseRing.N
	B := fastmath.NewRandomPolyMatrix(tau*d, l, ternarySampler, baseRing)
	BNTT := fastmath.NewPolyMatrix(tau*d, l, baseRing).NTT()
	BNTT.PopulateRows(func(i int) *fastmath.PolyNTTVec {
		nttRow := fastmath.NewPolyVec(l, baseRing).NTT()
		nttRow.Populate(func(j int) *fastmath.PolyNTT {
			return fastmath.ForceNTT(B.Get(i, j))
		})
		return nttRow
	})
	return BNTT
}

func CreateABPMask(tau int, ternarySampler fastmath.PolySampler, baseRing *ring.Ring) *fastmath.PolyNTTVec {
	abpMask := fastmath.NewRandomPolyVec(tau, ternarySampler, baseRing)
	abpMaskNTT := fastmath.NewPolyVec(tau, baseRing).NTT()
	for i := 0; i < tau; i++ {
		abpMaskNTT.Set(i, fastmath.ForceNTT(abpMask.Get(i)))
	}
	return abpMaskNTT
}

func CreateABPMaskedOpening(abpChal *fastmath.PolyNTTMatrix, abpMask *fastmath.PolyNTTVec, u *fastmath.IntVec, baseRing *ring.Ring) *fastmath.PolyNTTVec {
	d := baseRing.N
	// z = Bu + y
	uPadded := u
	if uPadded.Size() < abpChal.Cols()*d {
		uPadded = uPadded.Copy().Append(fastmath.NewIntVec(abpChal.Cols()*d-uPadded.Size(), baseRing))
	}
	B := abpChal.ToIntMatrix()
	z := B.MulVec(uPadded)
	z.Add(abpMask.ToIntVec())
	return z.UnderlyingPolysAsPolyNTTVec()
}

func NewABPEquation(abpChal *fastmath.PolyNTTMatrix, A *fastmath.IntMatrix, sIndex int, abpMask, abpMaskedOpening *fastmath.PolyNTTVec, baseRing *ring.Ring) *LinearEquation {
	fmt.Println("Creating ABP equation")
	d := baseRing.N
	// BAs + y = z
	APadded := A
	if APadded.Rows() < abpChal.Cols()*d {
		APadded = APadded.Copy().ExtendRows(fastmath.NewIntMatrix(abpChal.Cols()*d-APadded.Rows(), APadded.Cols(), baseRing))
	}
	BA := abpChal.ToIntMatrix().MulMat(APadded)
	y := abpMask.ToIntVec()
	z := abpMaskedOpening.ToIntVec()
	eqn := NewLinearEquation(z, BA.Cols())
	eqn.AppendDependentTerm(BA, sIndex)
	eqn.AppendVecTerm(y, baseRing)
	fmt.Println("Created ABP equation")
	return eqn
}
