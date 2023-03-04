package fastmath

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/ring"
)

// RingParams represents a polynomial ring.
type RingParams struct {
	BaseRing *ring.Ring
	Q        *big.Int
	D        int
	LogD     int
	Sigma    float64
	T        uint64
}

func BFVFullRingPN13() RingParams {
	// Initialize the ring parameters.
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParams, err := bfv.NewParametersFromLiteral(ringParamDef)
	if err != nil {
		panic("could not initialize the ring parameters: " + err.Error())
	}
	baseRing := ringParams.RingQP().RingQ
	return RingParams{
		BaseRing: baseRing,
		Q:        ringParams.Parameters.QBigInt(),
		D:        ringParams.N(),
		LogD:     ringParams.LogN(),
		Sigma:    ringParams.Sigma(),
		T:        ringParams.T(),
	}
}

func BFVFullRingCustom(paramLiteral bfv.ParametersLiteral, overrideLogD int) RingParams {
	// Initialize the ring parameters.
	ringParamDef := paramLiteral
	ringParamDef.LogN = overrideLogD
	ringParams, err := bfv.NewParametersFromLiteral(ringParamDef)
	if err != nil {
		panic("could not initialize the ring parameters: " + err.Error())
	}
	baseRing := ringParams.RingQP().RingQ
	return RingParams{
		BaseRing: baseRing,
		Q:        ringParams.Parameters.QBigInt(),
		D:        ringParams.N(),
		LogD:     ringParams.LogN(),
		Sigma:    ringParams.Sigma(),
		T:        ringParams.T(),
	}
}

func BFVFullShortCommtRing(logD int) RingParams {
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParamDef.LogN = logD
	ringParams, err := bfv.NewParametersFromLiteral(ringParamDef)
	if err != nil {
		panic("could not initialize the ring parameters: " + err.Error())
	}
	baseRing := ringParams.RingQP().RingQ
	return RingParams{
		BaseRing: baseRing,
		Q:        ringParams.QBigInt(),
		D:        ringParams.N(),
		LogD:     ringParams.LogN(),
		Sigma:    ringParams.Sigma(),
		T:        ringParams.T(),
	}
}

func BFVTwoLevelRingCustom(paramLiteral bfv.ParametersLiteral, overrideLogD int) RingParams {
	// Initialize the ring parameters.
	ringParamDef := paramLiteral
	ringParamDef.LogN = overrideLogD
	ringParams, err := bfv.NewParametersFromLiteral(ringParamDef)
	if err != nil {
		panic("could not initialize the ring parameters: " + err.Error())
	}
	baseRing := ringParams.RingQP().RingQ
	baseRing.Modulus = baseRing.Modulus[0:2]
	baseRing.ModulusAtLevel = baseRing.ModulusAtLevel[0:2]
	return RingParams{
		BaseRing: baseRing,
		Q:        baseRing.ModulusAtLevel[1],
		D:        ringParams.N(),
		LogD:     ringParams.LogN(),
		Sigma:    ringParams.Sigma(),
		T:        ringParams.T(),
	}
}

func BFVZeroLevelRingPN13() RingParams {
	// Initialize the ring parameters.
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParams, err := bfv.NewParametersFromLiteral(ringParamDef)
	if err != nil {
		panic("could not initialize the ring parameters: " + err.Error())
	}
	baseRing := ringParams.RingQP().RingQ
	baseRing.Modulus = baseRing.Modulus[0:1]
	baseRing.ModulusAtLevel = baseRing.ModulusAtLevel[0:1]
	return RingParams{
		BaseRing: baseRing,
		Q:        big.NewInt(int64(baseRing.Modulus[0])),
		D:        ringParams.N(),
		LogD:     ringParams.LogN(),
		Sigma:    ringParams.Sigma(),
		T:        ringParams.T(),
	}
}

func BFVZeroLevelRingCustom(paramLiteral bfv.ParametersLiteral, overrideLogD int) RingParams {
	// Initialize the ring parameters.
	ringParamDef := paramLiteral
	ringParamDef.LogN = overrideLogD
	ringParams, err := bfv.NewParametersFromLiteral(ringParamDef)
	if err != nil {
		panic("could not initialize the ring parameters: " + err.Error())
	}
	baseRing := ringParams.RingQP().RingQ
	baseRing.Modulus = baseRing.Modulus[0:1]
	baseRing.ModulusAtLevel = baseRing.ModulusAtLevel[0:1]
	return RingParams{
		BaseRing: baseRing,
		Q:        big.NewInt(int64(baseRing.Modulus[0])),
		D:        ringParams.N(),
		LogD:     ringParams.LogN(),
		Sigma:    ringParams.Sigma(),
		T:        ringParams.T(),
	}
}

func BFVZeroLevelShortCommtRing(logD int) RingParams {
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParamDef.LogN = logD
	ringParams, err := bfv.NewParametersFromLiteral(ringParamDef)
	if err != nil {
		panic("could not initialize the ring parameters: " + err.Error())
	}
	baseRing := ringParams.RingQP().RingQ
	baseRing.Modulus = baseRing.Modulus[0:1]
	baseRing.ModulusAtLevel = baseRing.ModulusAtLevel[0:1]
	return RingParams{
		BaseRing: baseRing,
		Q:        big.NewInt(int64(baseRing.Modulus[0])),
		D:        ringParams.N(),
		LogD:     ringParams.LogN(),
		Sigma:    ringParams.Sigma(),
		T:        ringParams.T(),
	}
}
