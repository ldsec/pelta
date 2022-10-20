package ens20

import (
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/tuneinsight/lattigo/v4/ring"
	"math/big"
)

type Settings struct {
	d               int      // deg(X^d + 1), a power of two
	q               *big.Int // Rational prime mod
	m               int      // # rows
	n               int      // # cols
	k               int      // Repetition rate
	delta1          int      // Width of the uniform distribution
	lambda          int      // M-LWE dimension
	kappa           int      // M-SIS dimension
	beta            float64  // Norm limit
	numSplits       int      // n / d
	baseRing        *ring.Ring
	uniformSampler  PolySampler
	ternarySampler  PolySampler
	gaussianSampler PolySampler
}

type PublicParams struct {
	A  math.Matrix
	u  math.Vector
	B0 math.Matrix
	b  math.Matrix
}

type Commitments struct {
	t0 math.Vector
}
