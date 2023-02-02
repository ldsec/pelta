package fastens20

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

// ProtocolConfig contains the settings for the augmented ENS20 protocol.
type ProtocolConfig struct {
	fastmath.RingParams
	M            int
	N            int
	Delta1       uint64
	K            int            // replication degree
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
	Cache *Cache
}

func NewProtocolConfig(ringParams fastmath.RingParams, delta1 uint64, numRows, numCols int) ProtocolConfig {
	uni, ter, gau := fastmath.GetSamplers(ringParams, delta1)
	return ProtocolConfig{
		RingParams:      ringParams,
		M:               numRows,
		N:               numCols,
		Delta1:          delta1,
		K:               1,
		Lambda:          1,
		Kappa:           1,
		TernarySlice:    fastmath.NewSlice(0, numCols),
		ABPEnabled:      false,
		Tau:             128,
		Bound:           big.NewInt(0),
		BoundSlice:      fastmath.NewSlice(0, 0),
		UniformSampler:  uni,
		TernarySampler:  ter,
		GaussianSampler: gau,
		Cache:           NewEmptyCache(),
	}
}

func (c ProtocolConfig) WithTernarySlice(ternarySlice fastmath.Slice) ProtocolConfig {
	c.TernarySlice = ternarySlice
	return c
}

func (c ProtocolConfig) WithReplication(k int) ProtocolConfig {
	c.K = k
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
	// TODO: fix
	return int(c.Delta1)
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

func GeneratePublicParameters(config ProtocolConfig, targetRel *crypto.ImmutLinearRelation) PublicParams {
	bSize := config.NumSplits() + 3
	B0 := fastmath.NewRandomPolyMatrix(config.Kappa,
		config.Lambda+config.Kappa+bSize,
		config.UniformSampler,
		config.BaseRing)
	b := fastmath.NewRandomPolyMatrix(bSize, B0.Cols(), config.UniformSampler, config.BaseRing)
	sig := fastmath.NewAutomorphism(uint64(config.D), uint64(config.K))
	// Build the cache.
	config.Cache.Build(config.K, sig, config.Q)
	return PublicParams{
		config: config,
		A:      targetRel.A,
		At:     targetRel.A.Transposed(),
		U:      targetRel.U,
		B0:     B0.NTT(),
		B:      b.NTT(),
		Sig:    sig,
	}
}
