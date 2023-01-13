package crypto

import (
	"fmt"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
	"github.com/tuneinsight/lattigo/v4/ring"
)

// CreateABPChallenge returns a random ternary polynomial matrix with dimensions tau x n.
func CreateABPChallenge(tau, n int, ternarySampler fastmath.PolySampler, baseRing *ring.Ring) *fastmath.IntMatrix {
	return fastmath.NewRandomIntMatrixFast(tau, n, ternarySampler, baseRing)
}

// CreateABPMask returns a random ternary vector of size tau.
// Warning: Returned y is in the NTT domain.
func CreateABPMask(tau int, ternarySampler fastmath.PolySampler, baseRing *ring.Ring) *fastmath.IntVec {
	// Create y of tau
	return fastmath.NewRandomIntVecFast(tau, ternarySampler, baseRing)
}

// CreateABPMaskedOpening returns R (abpChal) * s + y (abpMask)
func CreateABPMaskedOpening(abpChal *fastmath.IntMatrix, abpMask *fastmath.IntVec, s *fastmath.IntVec, baseRing *ring.Ring) *fastmath.IntVec {
	// z = Rs + y
	return abpChal.MulVec(s).Add(abpMask)
}

// NewABPEquation creates an equation of form Rs + y = z with s dependent where abpChal: R, abpMask: y, abpMaskedOpening: z.
func NewABPEquation(abpChal *fastmath.IntMatrix, sIndex int, abpMask, abpMaskedOpening *fastmath.IntVec, baseRing *ring.Ring) *LinearEquation {
	e := logging.LogExecStart("NewABPEquation", "creating")
	// Add the equation Rs + y = z
	R := abpChal.Copy()
	y := abpMask
	z := abpMaskedOpening
	if y.Size()%baseRing.N != 0 {
		padLength := baseRing.N - (y.Size() % baseRing.N)
		logging.Log("NewABPEquation",
			fmt.Sprintf("y (mask) is not a multiple of poly degree (size = %d), padding accordingly (+%d)", abpMask.Size(), padLength))
		y.Append(fastmath.NewIntVec(padLength, baseRing))
		z.Append(fastmath.NewIntVec(padLength, baseRing))
		R.ExtendRows(fastmath.NewIntMatrix(padLength, R.Cols(), baseRing))
	}
	eqn := NewLinearEquation(z, R.Cols())
	eqn.AppendDependentTerm(R, sIndex)
	eqn.AppendVecTerm(y, baseRing)
	e.LogExecEnd()
	return eqn
}
