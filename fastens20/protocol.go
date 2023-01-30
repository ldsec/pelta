package fastens20

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// ProtocolConfig contains the settings for the ENS20 protocol.
type ProtocolConfig struct {
	RingParams   fastmath.RingParams
	TargetRel    *crypto.ImmutLinearRelation
	BaseRing     *ring.Ring
	D            int      // poly. degree
	Q            *big.Int // ring modulus
	M            int      // # rows
	N            int      // # cols
	K            int      // replication degree
	InvK         uint64   // k^{-1} mod q
	Delta1       int
	Lambda       int            // security parameter
	Kappa        int            // security parameter
	TernarySlice fastmath.Slice // slice of s that should be ternary
	// abp
	ABPEnabled bool
	BoundSlice fastmath.Slice // slice of s that should be abp checked
	Tau        int            // abp security parameter
	Bound      *big.Int       // abp bound
	// samplers
	UniformSampler  fastmath.PolySampler
	TernarySampler  fastmath.PolySampler
	GaussianSampler fastmath.PolySampler
	// Cache
	ValueCache *Cache
}

// DefaultProtocolConfig returns the default configuration for an ENS20 execution.
// `ringParams` denotes the ring over which `relation` is defined.
func DefaultProtocolConfig(ringParams fastmath.RingParams, rel *crypto.ImmutLinearRelation) ProtocolConfig {
	// Create the samplers.
	prng, err := utils.NewPRNG()
	if err != nil {
		panic("could not initialize the prng")
	}
	// Construct the augmented ternary sampler.
	originalTernarySampler := ring.NewTernarySampler(prng, ringParams.BaseRing, 1.0/3.0, false)
	ternarySampler := fastmath.NewAugmentedTernarySampler(originalTernarySampler, ringParams.BaseRing)
	delta1 := 16
	numRows := rel.A.Rows()
	numCols := rel.A.Cols()
	invK := big.NewInt(0).ModInverse(big.NewInt(int64(1)), ringParams.Q).Uint64()
	return ProtocolConfig{
		RingParams:      ringParams,
		TargetRel:       rel,
		D:               ringParams.D,
		Q:               ringParams.Q,
		N:               numCols,
		M:               numRows,
		K:               1,
		InvK:            invK,
		Delta1:          delta1,
		Lambda:          1,
		Kappa:           1,
		TernarySlice:    fastmath.NewSlice(0, numCols),
		BaseRing:        ringParams.BaseRing,
		UniformSampler:  ring.NewUniformSampler(prng, ringParams.BaseRing),
		TernarySampler:  ternarySampler,
		GaussianSampler: ring.NewGaussianSampler(prng, ringParams.BaseRing, ringParams.Sigma, delta1),
		Tau:             128,
		Bound:           big.NewInt(0),
		ABPEnabled:      false,
		BoundSlice:      fastmath.NewSlice(0, 0),
		ValueCache:      NewEmptyCache(),
	}
}

func (c ProtocolConfig) WithTernarySlice(ternarySlice fastmath.Slice) ProtocolConfig {
	c.TernarySlice = ternarySlice
	return c
}

func (c ProtocolConfig) WithReplication(k int) ProtocolConfig {
	c.K = k
	c.InvK = big.NewInt(0).ModInverse(big.NewInt(int64(k)), c.Q).Uint64()
	return c
}

func (c ProtocolConfig) WithSecurityParameters(kappa, lambda int) ProtocolConfig {
	c.Kappa = kappa
	c.Lambda = lambda
	return c
}

func (c ProtocolConfig) WithABP(tau int, bound *big.Int, boundSlice fastmath.Slice) ProtocolConfig {
	c.ABPEnabled = true
	c.Tau = tau
	c.Bound = bound
	c.BoundSlice = boundSlice
	return c
}

// NumSplits returns n/d
func (c ProtocolConfig) NumSplits() int {
	return c.N / c.D
}

// NumTernarySplits returns the number of splits that should be checked to be ternary.
func (c ProtocolConfig) NumTernarySplits() int {
	return c.TernarySlice.Size() / c.D
}

// Beta returns the norm Gaussian limit.
func (c ProtocolConfig) Beta() int {
	return c.Delta1
}

// PublicParams contains the public parameters of the protocol.
type PublicParams struct {
	config ProtocolConfig
	A      fastmath.ImmutIntMatrix
	At     fastmath.ImmutIntMatrix
	U      *fastmath.IntVec
	B0     *fastmath.PolyNTTMatrix
	B      *fastmath.PolyNTTMatrix
	Sig    fastmath.Automorphism
}

func GeneratePublicParameters(config ProtocolConfig) PublicParams {
	bSize := config.NumSplits() + 3
	B0 := fastmath.NewRandomPolyMatrix(config.Kappa,
		config.Lambda+config.Kappa+bSize,
		config.UniformSampler,
		config.BaseRing)
	b := fastmath.NewRandomPolyMatrix(bSize, B0.Cols(), config.UniformSampler, config.BaseRing)
	sig := fastmath.NewAutomorphism(uint64(config.D), uint64(config.K))
	return PublicParams{
		config: config,
		A:      config.TargetRel.A,
		At:     config.TargetRel.A.Transposed(),
		U:      config.TargetRel.U,
		B0:     B0.NTT(),
		B:      b.NTT(),
		Sig:    sig,
	}
}