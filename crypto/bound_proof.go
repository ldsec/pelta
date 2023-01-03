package crypto

import (
	"fmt"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
)

// CreateABPChallenge returns a random ternary polynomial matrix with dimensions tau x m.
func CreateABPChallenge(tau, m int, ternarySampler fastmath.PolySampler, baseRing *ring.Ring) *fastmath.IntMatrix {
	// Create R of tau x m
	return fastmath.NewRandomTernaryIntMatrix(tau, m, baseRing)
}

// CreateABPMask returns a random ternary vector of size tau.
// Warning: Returned y is in the NTT domain.
func CreateABPMask(tau int, ternarySampler fastmath.PolySampler, baseRing *ring.Ring) *fastmath.IntVec {
	// Create y of tau
	return fastmath.NewRandomTernaryIntVec(tau, baseRing)
}

// CreateABPMaskedOpening returns abpChal * u + abpMask
func CreateABPMaskedOpening(abpChal *fastmath.IntMatrix, abpMask *fastmath.IntVec, u *fastmath.IntVec, baseRing *ring.Ring) *fastmath.IntVec {
	// z = Ru + y
	return abpChal.MulVec(u).Add(abpMask)
}

// NewABPEquation creates an equation of form RAs + y = z with s dependent where abpChal: R, abpMask: y, abpMaskedOpening: z.
func NewABPEquation(abpChal *fastmath.IntMatrix, A *fastmath.IntMatrix, sIndex int, abpMask, abpMaskedOpening *fastmath.IntVec, baseRing *ring.Ring) *LinearEquation {
	fmt.Println("creating ABP equation")
	// Add the equation RAs + y = z
	RA := abpChal.MulMat(A)
	y := abpMask
	z := abpMaskedOpening
	if y.Size()%baseRing.N != 0 {
		padLength := baseRing.N - (y.Size() % baseRing.N)
		fmt.Printf("y (mask) is not a multiple of poly degree (size = %d), padding accordingly (+%d)\n", abpMask.Size(), padLength)
		y.Append(fastmath.NewIntVec(padLength, baseRing))
		z.Append(fastmath.NewIntVec(padLength, baseRing))
		RA.ExtendRows(fastmath.NewIntMatrix(padLength, RA.Cols(), baseRing))
	}
	eqn := NewLinearEquation(z, RA.Cols())
	eqn.AppendDependentTerm(RA, sIndex)
	eqn.AppendVecTerm(y, baseRing)
	fmt.Println("created ABP equation")
	return eqn
}
