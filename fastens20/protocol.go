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

// DefaultProtocolConfig returns the default configuration for an ENS20 execution.
// `ringParams` denotes the ring over which `relation` is defined.
func DefaultProtocolConfig(ringParams fastmath.RingParams, rel crypto.LinearRelation) ProtocolConfig {
	// Create the samplers.
	prng, err := utils.NewPRNG()
	if err != nil {
		panic("could not initialize the prng: %s")
	}
	// Construct the augmented ternary sampler.
	originalTernarySampler := ring.NewTernarySampler(prng, ringParams.BaseRing, 1.0/3.0, false)
	ternarySampler := fastmath.NewAugmentedTernarySampler(originalTernarySampler, ringParams.BaseRing)
	delta1 := 16
	numRows := rel.A.Rows()
	numCols := rel.A.Cols()
	return ProtocolConfig{
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
		TernarySampler:  ternarySampler,
		GaussianSampler: ring.NewGaussianSampler(prng, ringParams.BaseRing, ringParams.Sigma, delta1),
	}
}

func (c ProtocolConfig) WithTernaryPrefix(ternaryLength int) ProtocolConfig {
	c.TernaryLength = ternaryLength
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

// NumSplits returns n/d
func (c ProtocolConfig) NumSplits() int {
	return c.N / c.D
}

// NumTernarySplits returns the number of splits that are in ternary.
func (c ProtocolConfig) NumTernarySplits() int {
	return c.TernaryLength / c.D
}

// Beta returns the norm Gaussian limit.
func (c ProtocolConfig) Beta() int {
	return c.Delta1
}

// PublicParams contains the public parameters of the protocol.
type PublicParams struct {
	config ProtocolConfig
	A      *fastmath.IntMatrix
	U      *fastmath.IntVec
	B0     *fastmath.PolyNTTMatrix
	B      *fastmath.PolyNTTMatrix
	Sig    fastmath.Automorphism
}

func GeneratePublicParameters(rel crypto.LinearRelation, config ProtocolConfig) PublicParams {
	bSize := config.NumSplits() + 3
	B0 := fastmath.NewRandomPolyMatrix(config.Kappa,
		config.Lambda+config.Kappa+bSize,
		config.UniformSampler,
		config.BaseRing)
	b := fastmath.NewRandomPolyMatrix(bSize, B0.Cols(), config.UniformSampler, config.BaseRing)
	sig := fastmath.NewAutomorphism(uint64(config.D), uint64(config.K))
	return PublicParams{
		config: config,
		A:      rel.A,
		U:      rel.U,
		B0:     B0.NTT(),
		B:      b.NTT(),
		Sig:    sig,
	}
}

func Execute(s *fastmath.IntVec, params PublicParams) bool {
	prover := NewProver(params)
	verifier := NewVerifier(params)
	// Commit to the message.
	t0, t, w, ps := prover.CommitToMessage(s)
	alpha, gamma, vs := verifier.CreateMasks(t0, t, w)
	t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
	c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
	// Recreate the masked opening until it satisfies the shortness condition.
	z, ps, err := prover.MaskedOpening(c, ps)
	for err != nil {
		z, ps, err = prover.MaskedOpening(c, ps)
	}
	return verifier.Verify(z, vs)
}
