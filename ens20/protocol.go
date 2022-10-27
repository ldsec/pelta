package ens20

import (
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
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
	BaseRing        *ring.Ring
	UniformSampler  math.PolySampler
	TernarySampler  math.PolySampler
	GaussianSampler math.PolySampler
}

// NumSplits returns n/d
func (s Settings) NumSplits() int {
	return s.N / s.D
}

// PublicParams contains the public parameters of the protocol.
type PublicParams struct {
	A  algebra.Matrix
	U  rings.IntVector
	B0 algebra.Matrix
	B  algebra.Matrix
}

func NewDummyPublicParameters(s rings.IntVector, settings Settings) PublicParams {
	A := math.NewRandomIntegerMatrix(settings.M, settings.N, settings.Q)
	// As = U
	u := rings.NewIntVec(A.MulVec(s.AsVec()))
	bSize := settings.NumSplits() + 3
	B0 := math.NewRandomPolynomialMatrix(settings.Kappa, settings.Lambda+settings.Kappa+bSize, settings.BaseRing, settings.UniformSampler)
	b := math.NewRandomPolynomialMatrix(bSize, B0.Cols(), settings.BaseRing, settings.UniformSampler)

	return PublicParams{
		A:  A,
		U:  u,
		B0: B0,
		B:  b,
	}
}
