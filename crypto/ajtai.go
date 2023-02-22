package crypto

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
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

// GetAjtaiKappa returns k s.t. (As + Br) - ([As + Br] mod p) = kp
func GetAjtaiKappa(comP, comQ *fastmath.IntVec, params AjtaiConfig) *fastmath.IntVec {
	// comQ - comP = (As + Br) - [(As + Br) mod p]
	diff := comQ.Copy().Add(comP.Copy().Neg())
	// kappa := (comQ - comP) / p s.t. kappa * p is the difference
	kappa := fastmath.NewIntVec(diff.Size(), params.RingParams.BaseRing)
	kappa.PopulateCoeffs(func(i int) fastmath.Coeff {
		diffi := diff.GetCoeff(i)
		for lvl := range diffi {
			qlvl := big.NewInt(int64(params.RingParams.BaseRing.Modulus[lvl]))
			pInv := big.NewInt(0).ModInverse(params.P, qlvl)
			res := big.NewInt(0).Mul(big.NewInt(int64(diffi[lvl])), pInv)
			res.Mod(res, qlvl)
			diffi[lvl] = res.Uint64()
		}
		return diffi
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
