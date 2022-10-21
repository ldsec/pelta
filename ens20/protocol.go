package ens20

import (
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math/big"
)

// Settings represents the protocol settings.
type Settings struct {
	D               int      // deg(X^D + 1), a power of two
	Q               *big.Int // Rational prime mod
	M               int      // # rows
	N               int      // # cols
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
	u  math.IntVector
	B0 math.Matrix
	b  math.Matrix
}
