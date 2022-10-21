package tests

import (
	"github.com/ldsec/codeBase/commitment/ens20"
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
	"testing"
)

func getTestSettings() ens20.Settings {
	// Initialize the ring parameters.
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParams, err := bfv.NewParametersFromLiteral(ringParamDef)
	if err != nil {
		panic("could not initialize the ring parameters: %s")
	}
	settings := ens20.Settings{
		D:      ringParams.N(),
		Q:      ringParams.QBigInt(),
		N:      ringParams.N(),
		M:      1,
		K:      1,
		Delta1: 16,
		Lambda: 1,
		Kappa:  1,
		Beta:   1,
	}
	// Initialize the ring.
	settings.BaseRing = ringParams.RingQP().RingQ
	// Create the samplers.
	prng, err := utils.NewPRNG()
	if err != nil {
		panic("could not initialize the prng: %s")
	}
	settings.UniformSampler = ring.NewUniformSampler(prng, settings.BaseRing)
	settings.TernarySampler = ring.NewTernarySampler(prng, settings.BaseRing, 1.0/3.0, false)
	settings.GaussianSampler = ring.NewGaussianSampler(prng, settings.BaseRing, ringParams.Sigma(), settings.Delta1)
	// Inputs
	settings.NumSplits = settings.N / settings.D
	return settings
}

func TestSimple(tst *testing.T) {
	settings := getTestSettings()
	s := ens20.NewRandomIntegerVector(settings.N, settings.Q)
	params := ens20.NewDummyPublicParameters(s, settings)
	prover := ens20.NewProver(params, settings)
	verifier := ens20.NewVerifier(params, settings)

	t0, t, w, ps := prover.CommitToMessage(s)
	alpha, gamma, vs := verifier.CreateMasks(t0, t, w)
	t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
	c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
	// Recreate the masked opening until it satisfies the shortness condition.
	z, ps, err := prover.MaskedOpening(c, ps)
	for err != nil {
		z, ps, err = prover.MaskedOpening(c, ps)
	}
	if !verifier.Verify(z, vs) {
		tst.Errorf("TestSimple: Verification failed!")
	}
}
