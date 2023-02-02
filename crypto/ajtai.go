package crypto

import (
	"github.com/ldsec/codeBase/commitment/fastmath"
	"math/big"
)

type AjtaiConfig struct {
	fastmath.RingParams
	P *big.Int // commitment modulo
}

func NewAjtaiConfig(p *big.Int, ringParams fastmath.RingParams) AjtaiConfig {
	return AjtaiConfig{ringParams, p}
}

// GetAjtaiCommitments returns (As + Br) and (As + Br) mod p
func GetAjtaiCommitments(A, B *fastmath.IntMatrix, s, r *fastmath.IntVec, config AjtaiConfig) (*fastmath.IntVec, *fastmath.IntVec) {
	comQ := A.MulVec(s).Add(B.MulVec(r))
	comP := comQ.Copy().Reduce(config.P)
	return comQ, comP
}

// GetKappa returns k s.t. (As + Br) - ([As + Br] mod p) = kp
func GetAjtaiKappa(comP, comQ *fastmath.IntVec, params AjtaiConfig) *fastmath.IntVec {
	// comQ - comP = (As + Br) - [(As + Br) mod p]
	diff := comQ.Copy().Add(comP.Copy().Neg())
	// kappa := (comQ - comP) / p s.t. kappa * p is the difference between
	// pInv := big.NewInt(0).ModInverse(p, config.Q).Uint64()
	kappa := fastmath.NewIntVec(diff.Size(), params.RingParams.BaseRing)
	kappa.Populate(func(i int) uint64 {
		return diff.GetLevel(i, 0) / params.P.Uint64()
	})
	return kappa
}

// NewPaddedAjtaiEquation returns the equation comP = As + Br - kp with rows padded up to the size of s. This is done to be able to generate
// the Id * (-p) on the lhs.
func NewPaddedAjtaiEquation(comP *fastmath.IntVec, A, B *fastmath.IntMatrix, s, r, kappa *fastmath.IntVec, config AjtaiConfig) *LinearEquation {
	// Num cols
	d := s.Size()
	// Num rows
	l := comP.Size()
	padLength := d - l
	paddedA := A.Copy().(*fastmath.IntMatrix)
	paddedA.ExtendRows(fastmath.NewIntMatrix(padLength, d, config.RingParams.BaseRing))
	paddedB := B.Copy().(*fastmath.IntMatrix)
	paddedB.ExtendRows(fastmath.NewIntMatrix(padLength, d, config.RingParams.BaseRing))
	paddedKappa := kappa.Copy().Append(fastmath.NewIntVec(padLength, config.RingParams.BaseRing))
	paddedComP := comP.Copy().Append(fastmath.NewIntVec(padLength, config.RingParams.BaseRing))
	negP := fastmath.NewCoeffFromBigInt(config.P, config.RingParams.BaseRing.Modulus).Neg(config.RingParams.BaseRing.Modulus)
	eqn := NewLinearEquation(paddedComP, d)
	eqn.AppendTerm(paddedA, s).
		AppendTerm(paddedB, r).
		AppendTerm(fastmath.NewIdIntMatrix(d, config.BaseRing).ScaleCoeff(negP), paddedKappa)
	return eqn
}
