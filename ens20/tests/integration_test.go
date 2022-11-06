package tests

import (
	"github.com/ldsec/codeBase/commitment/ens20"
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
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
	return settings
}

func ExecuteAndTestCorrectness(outputPrefix string, tst *testing.T, s rings.ZIntVector, settings ens20.Settings, params ens20.PublicParams) {
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

func ExecuteAndTestSoundness(outputPrefix string, tst *testing.T, s rings.ZIntVector, settings ens20.Settings, params ens20.PublicParams) {
	perturbArray := func(v *algebra.MultiArray) {
		v.ForEach(func(el algebra.Element, _ []int) {
			el.(rings.Polynomial).RRot(2)
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
	s := math.NewRandomTernaryIntegerVector(settings.N, settings.Q)
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
	if !t.Element(settings.NumSplits()).Copy().Add(params.B.Row(settings.NumSplits()).Copy().AsVec().Dot(ps.R).Neg()).Eq(ps.G) {
		tst.Errorf("TestConsistency: CommitToMessage t[n/d] != b[n/d] * r + g")
	}
	alpha, gamma, vs := verifier.CreateMasks(t0, t, w)
	t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
	c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
	// Recreate the masked opening until it satisfies the shortness condition.
	z, ps, err := prover.MaskedOpening(c, ps)
	for err != nil {
		z, ps, err = prover.MaskedOpening(c, ps)
	}
	// Constructing f
	f := algebra.NewMatrixFromDimensions(settings.K, settings.NumSplits()).Populate(
		func(i int, j int) algebra.Element {
			// t[j]*sig^i(c)
			tmp := t.Element(j).Copy().Mul(vs.Sig.Permute(int64(i), c))
			// (b[j] * z[i]) - (t[j]*sig^i(c))
			return params.B.Row(j).Copy().AsVec().Dot(z.Row(i)).Sub(tmp)
		})
	f2 := params.B.Row(settings.NumSplits() + 1).Copy().AsVec().
		Dot(z.Row(0)).
		Sub(c.Copy().
			Mul(t.Element(settings.NumSplits() + 1)))
	f3 := params.B.Row(settings.NumSplits() + 2).Copy().AsVec().
		Dot(z.Row(0)).
		Sub(c.Copy().
			Mul(t.Element(settings.NumSplits() + 2)))
	// Check some equalities.
	check1 := f.All(
		func(el algebra.Element, i int, j int) bool {
			rhs := params.B.Row(j).Copy().AsVec().Dot(ps.Y.Row(i)).
				Sub(ps.Sig.Permute(int64(i), c).Mul(ps.SHat.Element(j)))
			return el.Eq(rhs)
		})
	if !check1 {
		tst.Errorf("TestConsistency: f equivalence test failed")
	}
	check2Sum := ens20.CommitmentSum(settings.K, settings.NumSplits(), rings.NewPolyVec(alpha), ps.Sig,
		func(i int, j int) rings.Polynomial {
			// (b[j] * y[i])^2
			tmp := params.B.Row(j).Copy().AsVec().
				Dot(ps.Y.Row(i)).
				Pow(2)
			// 3s[j] - 3
			tmp2 := ps.SHat.Element(j).Copy().(rings.Polynomial).
				Scale(3).
				Sub(rings.NewOnePolynomial(settings.BaseRing).
					Scale(3))
			// (3s[j]-3) (b[j] * y[i])^2
			return tmp2.
				Mul(tmp).(rings.Polynomial)
		})
	check2 := f2.Eq(params.B.Row(settings.NumSplits() + 1).Copy().AsVec().
		Dot(ps.Y.Row(0)).
		Sub(params.B.Row(settings.NumSplits() + 2).Copy().AsVec().
			Dot(ps.Y.Row(0)).Mul(c)).
		Add(check2Sum.Mul(c)))
	if !check2 {
		tst.Errorf("TestConsistency: f[n/d+2] equivalence test failed")
	}
	check3Sum := ens20.CommitmentSum(settings.K, settings.NumSplits(), rings.NewPolyVec(alpha), ps.Sig,
		func(i int, j int) rings.Polynomial {
			// b[j]*y[i]
			tmp := params.B.Row(j).Copy().AsVec().
				Dot(ps.Y.Row(i))
			// 2s[j]-1
			tmp1 := ps.SHat.Element(j).Copy().Scale(2).
				Sub(rings.NewOnePolynomial(settings.BaseRing))
			// s[j]-2
			tmp2 := ps.SHat.Element(j).Copy().
				Sub(rings.NewOnePolynomial(settings.BaseRing).
					Scale(2))
			// (s[j]-1)s[j]
			tmp3 := ps.SHat.Element(j).Copy().
				Sub(rings.NewOnePolynomial(settings.BaseRing)).
				Mul(ps.SHat.Element(j))
			// [(2s[j]-1) (2s[j]-2) + s[j](s[j]-1)](b[j]*y[i])
			return tmp.Mul(tmp1.Mul(tmp2).Add(tmp3)).(rings.Polynomial)
		})
	check3 := f3.Eq(
		params.B.Row(settings.NumSplits() + 2).Copy().AsVec().
			Dot(ps.Y.Row(0)).
			Sub(check3Sum.Mul(c)))
	if !check3 {
		tst.Errorf("TestConsistency: f[n/d+3] equivalence test failed")
	}

}

func TestSimple(tst *testing.T) {
	settings := getSimpleTestSettings()
	s := math.NewRandomTernaryIntegerVector(settings.N, settings.Q)
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
	s := math.NewRandomTernaryIntegerVector(settings.N, settings.Q)
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
	s := math.NewRandomTernaryIntegerVector(settings.N, settings.Q)
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
	settings.K = 4
	s := math.NewRandomTernaryIntegerVector(settings.N, settings.Q)
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
	for _, k := range []int{1, 2, 4} {
		settings.K = k
		for _, numSplits := range []int{1, 5, 10, 20, 30} {
			settings.N = settings.D * numSplits
			tst.Logf("Executing for k=%d, n/d=%d...", k, settings.NumSplits())
			s := math.NewRandomTernaryIntegerVector(settings.N, settings.Q)
			params := ens20.NewDummyPublicParameters(s, settings)
			t0 := time.Now()
			ExecuteAndTestCorrectness("Simple", tst, s, settings, params)
			dt := time.Since(t0).Milliseconds()
			tst.Logf(">> Took %d ms", dt)
		}
	}
}
