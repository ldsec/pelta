package crypto

import (
	"fmt"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
)

// CreateABPChallenge returns a random ternary polynomial matrix with dimensions tau x m.
func CreateABPChallenge(tau, m int, ternarySampler fastmath.PolySampler, baseRing *ring.Ring) *fastmath.PolyNTTMatrix {
	// Create R^T of m x (tau * d)
	RT := fastmath.NewRandomPolyMatrix(m, tau, ternarySampler, baseRing)
	RTNTT := fastmath.NewPolyMatrix(RT.Rows(), RT.Cols(), baseRing).NTT()
	RTNTT.PopulateRows(func(i int) *fastmath.PolyNTTVec {
		nttRow := fastmath.NewPolyVec(RT.Cols(), baseRing).NTT()
		nttRow.Populate(func(j int) *fastmath.PolyNTT {
			return fastmath.ForceNTT(RT.Get(i, j))
		})
		return nttRow
	})
	return RTNTT
}

// CreateABPMask returns a random ternary vector of size tau.
func CreateABPMask(tau int, ternarySampler fastmath.PolySampler, baseRing *ring.Ring) *fastmath.PolyNTTVec {
	// Create y of (tau * d)
	abpMask := fastmath.NewRandomPolyVec(tau, ternarySampler, baseRing)
	abpMaskNTT := fastmath.NewPolyVec(tau, baseRing).NTT()
	abpMaskNTT.Populate(func(i int) *fastmath.PolyNTT {
		return fastmath.ForceNTT(abpMask.Get(i))
	})
	return abpMaskNTT
}

// CreateABPMaskedOpening returns abpChal * u + abpMask
func CreateABPMaskedOpening(abpChal *fastmath.PolyNTTMatrix, abpMask *fastmath.PolyNTTVec, u *fastmath.IntVec, baseRing *ring.Ring) *fastmath.PolyNTTVec {
	// z = Ru + y
	R := abpChal.ToIntMatrix().Transposed()
	z := R.MulVec(u)
	z.Add(abpMask.ToIntVec())
	return z.UnderlyingPolysAsPolyNTTVec()
}

// NewABPEquation creates an equation of form RAs + y = z with s dependent where abpChal: R, abpMask: y, abpMaskedOpening: z.
func NewABPEquation(abpChal *fastmath.PolyNTTMatrix, A *fastmath.IntMatrix, sIndex int, abpMask, abpMaskedOpening *fastmath.PolyNTTVec, baseRing *ring.Ring) *LinearEquation {
	fmt.Println("creating ABP equation (1)")
	// Add the equation RAs + y = z
	RA := abpChal.ToIntMatrix().Transposed().MulMat(A)
	y := abpMask.ToIntVec()
	z := abpMaskedOpening.ToIntVec()
	fmt.Println("creating ABP equation (2)")
	eqn := NewLinearEquation(z, RA.Cols())
	eqn.AppendDependentTerm(RA, sIndex)
	eqn.AppendVecTerm(y, baseRing)
	fmt.Println("created ABP equation")
	return eqn
}
