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
	BaseRing        *ring.Ring
	UniformSampler  fastmath.PolySampler
	TernarySampler  fastmath.PolySampler
	GaussianSampler fastmath.PolySampler
}

// DefaultConfig returns the default configuration for an ENS20 execution.
func DefaultConfig(ringParams fastmath.RingParams, numRows int, numCols int) Config {
	// Create the samplers.
	prng, err := utils.NewPRNG()
	if err != nil {
		panic("could not initialize the prng: %s")
	}
	delta1 := 16
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
	Sig    fastmath.Automorphism
}

func GeneratePublicParameters(sis crypto.SISProblem, protocolConfig Config) PublicParams {
	bSize := protocolConfig.NumSplits() + 3
	B0 := fastmath.NewRandomPolyMatrix(protocolConfig.Kappa,
		protocolConfig.Lambda+protocolConfig.Kappa+bSize,
		protocolConfig.UniformSampler,
		protocolConfig.BaseRing)
	b := fastmath.NewRandomPolyMatrix(bSize, B0.Cols(), protocolConfig.UniformSampler, protocolConfig.BaseRing)
	sig := fastmath.NewAutomorphism(uint64(protocolConfig.D), uint64(protocolConfig.K))
	return PublicParams{
		config: protocolConfig,
		A:      sis.A,
		U:      sis.U,
		B0:     B0.NTT(),
		B:      b.NTT(),
		Sig:    sig,
	}
}
