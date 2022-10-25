package tests

import (
	"github.com/ldsec/codeBase/commitment/ens20"
	"github.com/ldsec/codeBase/commitment/math"
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
		Q:      ringParams.RingQP().RingQ.ModulusAtLevel[0],
		N:      ringParams.N(),
		M:      16,
		K:      1,
		Delta1: 16,
		Lambda: 1,
		Kappa:  1,
		Beta:   16,
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

func ExecuteAndTest(outputPrefix string, tst *testing.T, s math.IntVector, settings ens20.Settings, params ens20.PublicParams) {
	prover := ens20.NewProver(params, settings)
	verifier := ens20.NewVerifier(params, settings)
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
	if !verifier.Verify(z, vs) {
		tst.Errorf(outputPrefix + ".TestProtocol: Verification failed!")
	}
}

func TestSimple(tst *testing.T) {
	settings := getTestSettings()
	s := ens20.NewRandomTernaryIntegerVector(settings.N, settings.Q)
	//fmt.Println("s = " + s.String())
	params := ens20.NewDummyPublicParameters(s, settings)
	prover := ens20.NewProver(params, settings)
	verifier := ens20.NewVerifier(params, settings)
	// Commit to the message.
	t0, t, w, ps := prover.CommitToMessage(s)
	// Check consistency in the state.
	if !ps.T0.Eq(t0.MultiArray) || !ps.T.Eq(t.MultiArray) || !ps.W.Eq(w.MultiArray) {
		tst.Errorf("TestSimple: CommitToMessage state consistency check failed")
	}
	// Check t0.
	if !t0.Eq(params.B0.Copy().AsMatrix().MulVec(ps.R).MultiArray) {
		tst.Errorf("TestSimple: CommitToMessage t0 != B0 * r")
	}
	// Check t[n/d].
	if !t.Element(settings.NumSplits).Copy().Add(params.B.Row(settings.NumSplits).Copy().AsVec().Dot(ps.R).Neg()).Eq(ps.G) {
		tst.Errorf("TestSimple: CommitToMessage t[n/d] != b[n/d] * r + g")
	}
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

func TestMultiReplication(tst *testing.T) {
	settings := getTestSettings()
	settings.K = 4
	s := ens20.NewRandomTernaryIntegerVector(settings.N, settings.Q)
	//fmt.Println("s = " + s.String())
	params := ens20.NewDummyPublicParameters(s, settings)
	ExecuteAndTest("MultiReplication", tst, s, settings, params)
}

func TestMultiSplit(tst *testing.T) {
	settings := getTestSettings()
	settings.N *= 4
	settings.NumSplits = 4
	s := ens20.NewRandomTernaryIntegerVector(settings.N, settings.Q)
	//fmt.Println("s = " + s.String())
	params := ens20.NewDummyPublicParameters(s, settings)
	ExecuteAndTest("MultiReplication", tst, s, settings, params)
}

func TestMultiSplitMultiReplication(tst *testing.T) {
	settings := getTestSettings()
	settings.N *= 4
	settings.NumSplits = 4
	settings.K = 4
	s := ens20.NewRandomTernaryIntegerVector(settings.N, settings.Q)
	//fmt.Println("s = " + s.String())
	params := ens20.NewDummyPublicParameters(s, settings)
	ExecuteAndTest("MultiReplication", tst, s, settings, params)
}
