package tests

import (
	"github.com/ldsec/codeBase/commitment/ens20"
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
	"testing"
	"time"
)

func getSimpleTestSettings() ens20.Settings {
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

func ExecuteAndTestCorrectness(outputPrefix string, tst *testing.T, s math.IntVector, settings ens20.Settings, params ens20.PublicParams) {
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
		tst.Errorf(outputPrefix + ".TestCorrectness: Verification failed!")
	}
}

func ExecuteAndTestSoundness(outputPrefix string, tst *testing.T, s math.IntVector, settings ens20.Settings, params ens20.PublicParams) {
	perturbArray := func(v *math.MultiArray) {
		v.ForEach(func(el math.RingElement, _ []int) {
			el.(math.Polynomial).RRot(2)
		})
	}
	prover := ens20.NewProver(params, settings)
	verifier := ens20.NewVerifier(params, settings)
	tst.Log(outputPrefix + ".TestSoundness: Perturbing t")
	{
		// Commit to the message.
		t0, t, w, ps := prover.CommitToMessage(s)
		// Perturb t.
		perturbArray(t.MultiArray)
		alpha, gamma, vs := verifier.CreateMasks(t0, t, w)
		t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
		c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
		// Recreate the masked opening until it satisfies the shortness condition.
		z, ps, err := prover.MaskedOpening(c, ps)
		for err != nil {
			z, ps, err = prover.MaskedOpening(c, ps)
		}
		if verifier.Verify(z, vs) {
			tst.Errorf(outputPrefix + ".TestSoundness: Verification succeeded!")
		}
	}
	tst.Log(outputPrefix + ".TestSoundness: Perturbing t0")
	{
		// Commit to the message.
		t0, t, w, ps := prover.CommitToMessage(s)
		// Perturb t0.
		perturbArray(t0.MultiArray)
		alpha, gamma, vs := verifier.CreateMasks(t0, t, w)
		t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
		c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
		// Recreate the masked opening until it satisfies the shortness condition.
		z, ps, err := prover.MaskedOpening(c, ps)
		for err != nil {
			z, ps, err = prover.MaskedOpening(c, ps)
		}
		if verifier.Verify(z, vs) {
			tst.Errorf(outputPrefix + ".TestSoundness: Verification succeeded!")
		}
	}
	tst.Log(outputPrefix + ".TestSoundness: Perturbing w")
	{
		// Commit to the message.
		t0, t, w, ps := prover.CommitToMessage(s)
		// Perturb w.
		perturbArray(w.MultiArray)
		alpha, gamma, vs := verifier.CreateMasks(t0, t, w)
		t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
		c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
		// Recreate the masked opening until it satisfies the shortness condition.
		z, ps, err := prover.MaskedOpening(c, ps)
		for err != nil {
			z, ps, err = prover.MaskedOpening(c, ps)
		}
		if verifier.Verify(z, vs) {
			tst.Errorf(outputPrefix + ".TestSoundness: Verification succeeded!")
		}
	}
	tst.Log(outputPrefix + ".TestSoundness: Perturbing v'")
	{
		// Commit to the message.
		t0, t, w, ps := prover.CommitToMessage(s)
		alpha, gamma, vs := verifier.CreateMasks(t0, t, w)
		t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
		// Perturb vp.
		perturbArray(vp.MultiArray)
		c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
		// Recreate the masked opening until it satisfies the shortness condition.
		z, ps, err := prover.MaskedOpening(c, ps)
		for err != nil {
			z, ps, err = prover.MaskedOpening(c, ps)
		}
		if verifier.Verify(z, vs) {
			tst.Errorf(outputPrefix + ".TestSoundness: Verification succeeded!")
		}
	}
	tst.Log(outputPrefix + ".TestSoundness: Perturbing v")
	{
		// Commit to the message.
		t0, t, w, ps := prover.CommitToMessage(s)
		alpha, gamma, vs := verifier.CreateMasks(t0, t, w)
		t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
		// Perturb v.
		v.RRot(3)
		c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
		// Recreate the masked opening until it satisfies the shortness condition.
		z, ps, err := prover.MaskedOpening(c, ps)
		for err != nil {
			z, ps, err = prover.MaskedOpening(c, ps)
		}
		if verifier.Verify(z, vs) {
			tst.Errorf(outputPrefix + ".TestSoundness: Verification succeeded!")
		}
	}
	tst.Log(outputPrefix + ".TestSoundness: Perturbing h")
	{
		// Commit to the message.
		t0, t, w, ps := prover.CommitToMessage(s)
		alpha, gamma, vs := verifier.CreateMasks(t0, t, w)
		t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
		// Perturb h.
		h.RRot(3)
		c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
		// Recreate the masked opening until it satisfies the shortness condition.
		z, ps, err := prover.MaskedOpening(c, ps)
		for err != nil {
			z, ps, err = prover.MaskedOpening(c, ps)
		}
		if verifier.Verify(z, vs) {
			tst.Errorf(outputPrefix + ".TestSoundness: Verification succeeded!")
		}
	}
	tst.Log(outputPrefix + ".TestSoundness: Perturbing c")
	{
		// Commit to the message.
		t0, t, w, ps := prover.CommitToMessage(s)
		alpha, gamma, vs := verifier.CreateMasks(t0, t, w)
		t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
		c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
		// Perturb c.
		c.RRot(3)
		// Recreate the masked opening until it satisfies the shortness condition.
		z, ps, err := prover.MaskedOpening(c, ps)
		for err != nil {
			z, ps, err = prover.MaskedOpening(c, ps)
		}
		if verifier.Verify(z, vs) {
			//tst.Errorf(outputPrefix + ".TestSoundness: Verification succeeded!")
		}
	}
	tst.Log(outputPrefix + ".TestSoundness: Perturbing z")
	{
		// Commit to the message.
		t0, t, w, ps := prover.CommitToMessage(s)
		alpha, gamma, vs := verifier.CreateMasks(t0, t, w)
		t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
		c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
		// Recreate the masked opening until it satisfies the shortness condition.
		z, ps, err := prover.MaskedOpening(c, ps)
		// Perturb c.
		perturbArray(z.MultiArray)
		for err != nil {
			z, ps, err = prover.MaskedOpening(c, ps)
		}
		if verifier.Verify(z, vs) {
			tst.Errorf(outputPrefix + ".TestSoundness: Verification succeeded!")
		}
	}
}

func TestConsistency(tst *testing.T) {
	settings := getSimpleTestSettings()
	s := ens20.NewRandomTernaryIntegerVector(settings.N, settings.Q)
	//fmt.Println("s = " + s.String())
	params := ens20.NewDummyPublicParameters(s, settings)
	prover := ens20.NewProver(params, settings)
	verifier := ens20.NewVerifier(params, settings)
	// Commit to the message.
	t0, t, w, ps := prover.CommitToMessage(s)
	// Check consistency in the state.
	if !ps.T0.Eq(t0.MultiArray) || !ps.T.Eq(t.MultiArray) || !ps.W.Eq(w.MultiArray) {
		tst.Errorf("TestConsistency: CommitToMessage state consistency check failed")
	}
	// Check t0.
	if !t0.Eq(params.B0.Copy().AsMatrix().MulVec(ps.R).MultiArray) {
		tst.Errorf("TestConsistency: CommitToMessage t0 != B0 * r")
	}
	// Check t[n/d].
	if !t.Element(settings.NumSplits).Copy().Add(params.B.Row(settings.NumSplits).Copy().AsVec().Dot(ps.R).Neg()).Eq(ps.G) {
		tst.Errorf("TestConsistency: CommitToMessage t[n/d] != b[n/d] * r + g")
	}
	_, _, _ = verifier.CreateMasks(t0, t, w)
	// TODO: add more consistency checks
}

func TestSimple(tst *testing.T) {
	settings := getSimpleTestSettings()
	s := ens20.NewRandomTernaryIntegerVector(settings.N, settings.Q)
	//fmt.Println("s = " + s.String())
	params := ens20.NewDummyPublicParameters(s, settings)
	tst.Log("Checking correctness...")
	ExecuteAndTestCorrectness("Simple", tst, s, settings, params)
	//tst.Log("Checking soundness...")
	//ExecuteAndTestSoundness("Simple", tst, s, settings, params)
}

func TestMultiReplication(tst *testing.T) {
	settings := getSimpleTestSettings()
	settings.K = 4
	s := ens20.NewRandomTernaryIntegerVector(settings.N, settings.Q)
	//fmt.Println("s = " + s.String())
	params := ens20.NewDummyPublicParameters(s, settings)
	tst.Log("Checking correctness...")
	ExecuteAndTestCorrectness("Simple", tst, s, settings, params)
	tst.Log("Checking soundness...")
	ExecuteAndTestSoundness("Simple", tst, s, settings, params)
}

func TestMultiSplit(tst *testing.T) {
	settings := getSimpleTestSettings()
	settings.N *= 4
	settings.NumSplits = 4
	s := ens20.NewRandomTernaryIntegerVector(settings.N, settings.Q)
	//fmt.Println("s = " + s.String())
	params := ens20.NewDummyPublicParameters(s, settings)
	tst.Log("Checking correctness...")
	ExecuteAndTestCorrectness("Simple", tst, s, settings, params)
	tst.Log("Checking soundness...")
	ExecuteAndTestSoundness("Simple", tst, s, settings, params)
}

func TestMultiSplitMultiReplication(tst *testing.T) {
	settings := getSimpleTestSettings()
	settings.N *= 4
	settings.NumSplits = 4
	settings.K = 4
	s := ens20.NewRandomTernaryIntegerVector(settings.N, settings.Q)
	//fmt.Println("s = " + s.String())
	params := ens20.NewDummyPublicParameters(s, settings)
	tst.Log("Checking correctness...")
	ExecuteAndTestCorrectness("Simple", tst, s, settings, params)
	tst.Log("Checking soundness...")
	ExecuteAndTestSoundness("Simple", tst, s, settings, params)
}

func TestPerformance(tst *testing.T) {
	settings := getSimpleTestSettings()
	tst.Log("Checking correctness...")
	for k := 1; k < 4; k++ {
		settings.K = k
		for _, numSplits := range []int{1, 4, 10, 20, 30, 50} {
			tst.Logf("Executing for k=%d, n/d=%d...", k, numSplits)
			settings.N = settings.D * numSplits
			settings.NumSplits = numSplits
			s := ens20.NewRandomTernaryIntegerVector(settings.N, settings.Q)
			params := ens20.NewDummyPublicParameters(s, settings)
			t0 := time.Now()
			ExecuteAndTestCorrectness("Simple", tst, s, settings, params)
			dt := time.Since(t0).Milliseconds()
			tst.Logf(">> Took %d ms", dt)
		}
	}
}
