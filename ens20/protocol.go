package ens20

import (
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math/big"
)

// Settings represents the protocol settings.
type Settings struct {
	D               int      // deg(X^D + 1), a power of two
	Q               *big.Int // Prime mod
	M               int      // # rows
	N               int      // # cols, must be >= D
	K               int      // Repetition rate
	Delta1          int      // Width of the uniform distribution
	Lambda          int      // M-LWE dimension
	Kappa           int      // M-SIS dimension
	Beta            float64  // Norm limit
	NumSplits       int      // N / D
	BaseRing        *ring.Ring
	UniformSampler  PolySampler
	TernarySampler  PolySampler
	GaussianSampler PolySampler
}

// PublicParams contains the public parameters of the protocol.
type PublicParams struct {
	A  math.Matrix
	U  math.IntVector
	B0 math.Matrix
	B  math.Matrix
}

func NewDummyPublicParameters(s math.IntVector, settings Settings) PublicParams {
	A := NewRandomIntegerMatrix(settings.M, settings.N, settings.Q)
	// As = U
	u := A.MulVec(s.AsVec()).AsIntVec()
	bSize := settings.NumSplits + 3
	B0 := NewRandomPolynomialMatrix(settings.Kappa, settings.Lambda+settings.Kappa+bSize, settings.BaseRing, settings.UniformSampler)
	b := NewRandomPolynomialMatrix(bSize, B0.Cols(), settings.BaseRing, settings.UniformSampler)

	return PublicParams{
		A:  A,
		U:  u,
		B0: B0,
		B:  b,
	}
}
