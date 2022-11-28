package fastens20

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Config contains the settings for the ENS20 protocol.
type Config struct {
	RingParams      fastmath.RingParams
	D               int
	Q               *big.Int
	M               int
	N               int
	K               int
	Delta1          int
	Lambda          int
	Kappa           int
	TernaryLength   int
	NormCheck       bool
	BaseRing        *ring.Ring
	UniformSampler  fastmath.PolySampler
	TernarySampler  fastmath.PolySampler
	GaussianSampler fastmath.PolySampler
}

// DefaultConfig returns the default configuration for an ENS20 execution.
// `ringParams` denotes the ring over which `relation` is defined.
func DefaultConfig(ringParams fastmath.RingParams, rel crypto.LinearRelation) Config {
	// Create the samplers.
	prng, err := utils.NewPRNG()
	if err != nil {
		panic("could not initialize the prng: %s")
	}
	delta1 := 16
	numRows := rel.A.Rows()
	numCols := rel.A.Cols()
	return Config{
		RingParams:      ringParams,
		D:               ringParams.D,
		Q:               ringParams.Q,
		N:               numCols,
		M:               numRows,
		K:               1,
		Delta1:          delta1,
		Lambda:          1,
		Kappa:           1,
		TernaryLength:   numCols,
		NormCheck:       false,
		BaseRing:        ringParams.BaseRing,
		UniformSampler:  ring.NewUniformSampler(prng, ringParams.BaseRing),
		TernarySampler:  ring.NewTernarySampler(prng, ringParams.BaseRing, 1.0/3.0, false),
		GaussianSampler: ring.NewGaussianSampler(prng, ringParams.BaseRing, ringParams.Sigma, delta1),
	}
}

func (c Config) WithTernaryPrefix(ternaryLength int) Config {
	c.TernaryLength = ternaryLength
	return c
}

func (c Config) WithReplication(k int) Config {
	c.K = k
	return c
}

func (c Config) WithSecurityParameters(kappa, lambda int) Config {
	c.Kappa = kappa
	c.Lambda = lambda
	return c
}

func (c Config) WithNormCheck() Config {
	c.NormCheck = true
	return c
}

// NumSplits returns n/d
func (c Config) NumSplits() int {
	return c.N / c.D
}

// NumTernarySplits returns the number of splits that are in ternary.
func (c Config) NumTernarySplits() int {
	return c.TernaryLength / c.D
}

// Beta returns the norm Gaussian limit.
func (c Config) Beta() int {
	return c.Delta1
}

// PublicParams contains the public parameters of the protocol.
type PublicParams struct {
	config Config
	A      *fastmath.IntMatrix
	U      *fastmath.IntVec
	B0     *fastmath.PolyNTTMatrix
	B      *fastmath.PolyNTTMatrix
	Bs     *fastmath.PolyNTTVec
	Sig    fastmath.Automorphism
}

func GeneratePublicParameters(rel crypto.LinearRelation, config Config) PublicParams {
	bSize := config.NumSplits() + 3
	B0 := fastmath.NewRandomPolyMatrix(config.Kappa,
		config.Lambda+config.Kappa+bSize,
		config.UniformSampler,
		config.BaseRing)
	b := fastmath.NewRandomPolyMatrix(bSize, B0.Cols(), config.UniformSampler, config.BaseRing)
	bs := fastmath.NewRandomPolyVec(B0.Cols(), config.UniformSampler, config.BaseRing)
	sig := fastmath.NewAutomorphism(uint64(config.D), uint64(config.K))
	return PublicParams{
		config: config,
		A:      rel.A,
		U:      rel.U,
		B0:     B0.NTT(),
		B:      b.NTT(),
		Bs:     bs.NTT(),
		Sig:    sig,
	}
}
