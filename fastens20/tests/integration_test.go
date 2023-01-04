package tests

import (
	"math/big"
	"testing"
	"time"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func ExecuteAndTestCorrectness(tst *testing.T, s *fastmath.IntVec, params fastens20.PublicParams) {
	prover := fastens20.NewProver(params)
	verifier := fastens20.NewVerifier(params)
	// Commit to the message.
	start_t := time.Now()
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
		tst.Errorf("verification failed")
	}
	end_t := time.Now()
	delta_t := end_t.Sub(start_t)
	tst.Logf("protocol execution took %dms", delta_t.Milliseconds())
}

func TestConsistency(tst *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 16
	n := bfvRing.D
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(n, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(rel, config)
	ACopy := params.A.Copy()
	BCopy := params.B.Copy()
	B0Copy := params.B0.Copy()
	prover := fastens20.NewProver(params)
	verifier := fastens20.NewVerifier(params)
	// Commit to the message.
	t0, t, w, ps := prover.CommitToMessage(s)
	// Check consistency in the state.
	if !ps.T0.Eq(t0) || !ps.T.Eq(t) || !ps.W.Eq(w) {
		tst.Errorf("CommitToMessage state consistency check failed")
	}
	// Check t0.
	if !t0.Eq(params.B0.Copy().MulVec(ps.R)) {
		tst.Errorf("CommitToMessage t0 != B0 * r")
	}
	// Check t[n/d].
	if !t.Get(config.NumSplits()).Copy().Add(params.B.Row(config.NumSplits()).Copy().Dot(ps.R).Neg()).Eq(ps.G) {
		tst.Errorf("CommitToMessage t[n/d] != b[n/d] * r + g")
	}
	alpha, gamma, vs := verifier.CreateMasks(t0, t, w)
	t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
	c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
	// Recreate the masked opening until it satisfies the shortness condition.
	z, ps, err := prover.MaskedOpening(c, ps)
	for err != nil {
		z, ps, err = prover.MaskedOpening(c, ps)
	}
	// End of the protocol.
	// Check that the parameters are not changed.
	if !ACopy.Eq(params.A) || !BCopy.Eq(params.B) || !B0Copy.Eq(params.B0) {
		tst.Errorf("parameters altered")
	}
	// Check that the states are consistent.
	if !ps.T.Eq(t) || !ps.T0.Eq(t0) || !ps.H.Eq(h) || !ps.V.Eq(v) || !ps.Vp.Eq(vp) || !ps.W.Eq(w) {
		tst.Errorf("prover state altered")
	}
	// Constructing f
	f := fastmath.NewPolyMatrix(config.K, config.NumSplits(), config.BaseRing).NTT()
	f.Populate(func(i, j int) *fastmath.PolyNTT {
		// t[j]*sig^i(c)
		tmp := t.Get(j).Copy().Mul(c.Permute(int64(i), params.Sig).NTT())
		// (b[j] * z[i]) - (t[j]*sig^i(c))
		return params.B.Row(j).Dot(z.Row(i)).Add(tmp.Neg())
	})
	// Check some equalities.
	check1 := f.All(func(i, j int, el *fastmath.PolyNTT) bool {
		// (b[j] * y[i]) - (s[j]*sig^i(c))
		rhs := params.B.Row(j).Copy().Dot(ps.Y.Row(i)).
			Add(c.Permute(int64(i), params.Sig).NTT().Mul(ps.SHat.Get(j)).Neg())
		return el.Eq(rhs)
	})
	if !check1 {
		tst.Errorf("f equivalence test failed")
	}
	// Check f[n/d + 2]
	f2 := params.B.Row(config.NumSplits() + 1).Copy().
		Dot(z.Row(0)).
		Add(c.Copy().NTT().
			Mul(t.Get(config.NumSplits() + 1)).Neg())
	check2Sum := fastens20.CommitmentSum(config.K, config.NumSplits(), alpha,
		func(i, j int) *fastmath.PolyNTT {
			// (b[j] * y[i])^2
			tmp := params.B.Row(j).Copy().
				Dot(ps.Y.Row(i)).
				Pow(2, config.Q.Uint64())
			// 3s[j] - 3
			tmp2 := ps.SHat.Get(j).Copy().
				Scale(3).
				Add(fastmath.NewOnePoly(3, config.BaseRing).Neg().NTT())
			// (3s[j]-3) (b[j] * y[i])^2
			return tmp2.Mul(tmp)
		}, params)
	check2RHS := params.B.Row(config.NumSplits() + 1).Copy().
		Dot(ps.Y.Row(0)).
		Add(params.B.Row(config.NumSplits() + 2).Copy().
			Dot(ps.Y.Row(0)).Mul(c.Copy().NTT()).Neg()).
		Add(check2Sum.Mul(c.Copy().NTT()))
	if !f2.Eq(check2RHS) {
		tst.Errorf("f[n/d+2] equivalence test failed")
	}
	// Check f[n/d + 3]
	f3 := params.B.Row(config.NumSplits() + 2).Copy().
		Dot(z.Row(0)).
		Add(c.Copy().NTT().
			Mul(t.Get(config.NumSplits() + 2)).Neg())
	check3Sum := fastens20.CommitmentSum(config.K, config.NumSplits(), alpha,
		func(i, j int) *fastmath.PolyNTT {
			// b[j]*y[i]
			tmp := params.B.Row(j).Dot(ps.Y.Row(i))
			// 2s[j]-1
			tmp1 := ps.SHat.Get(j).Copy().Scale(2).
				Add(fastmath.NewOnePoly(1, config.BaseRing).Neg().NTT())
			// s[j]-2
			tmp2 := ps.SHat.Get(j).Copy().
				Add(fastmath.NewOnePoly(2, config.BaseRing).Neg().NTT())
			// (s[j]-1)*s[j]
			tmp3 := ps.SHat.Get(j).Copy().
				Add(fastmath.NewOnePoly(1, config.BaseRing).Neg().NTT()).
				Mul(ps.SHat.Get(j))
			// [(2s[j]-1)*(s[j]-2) + s[j](s[j]-1)] * (b[j]*y[i])
			return tmp.Mul(tmp1.Mul(tmp2).Add(tmp3))
		}, params)
	check3RHS := params.B.Row(config.NumSplits() + 2).
		Dot(ps.Y.Row(0)).
		Add(check3Sum.Mul(c.Copy().NTT()).Neg())
	if !f3.Eq(check3RHS) {
		tst.Errorf("f[n/d+3] equivalence test failed")
	}

}

func TestSimple(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 16
	n := bfvRing.D
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(n, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(rel, config)
	ExecuteAndTestCorrectness(t, s, params)
}

func TestMultiReplication(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 16
	n := bfvRing.D
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(n, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel).
		WithReplication(4)
	params := fastens20.GeneratePublicParameters(rel, config)
	ExecuteAndTestCorrectness(t, s, params)
}

func TestMultiSplit(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 16
	n := bfvRing.D * 4
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(n, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(rel, config)
	ExecuteAndTestCorrectness(t, s, params)
}

func TestMultiSplitMultiReplication(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 16
	n := bfvRing.D * 4
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(n, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel).
		WithReplication(4)
	params := fastens20.GeneratePublicParameters(rel, config)
	ExecuteAndTestCorrectness(t, s, params)
}

func TestTernarySlice(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 16
	n := bfvRing.D * 4
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomIntVec(n, big.NewInt(5), bfvRing.BaseRing)
	ternarySlice := fastmath.NewSlice(n/2, n)
	for i := ternarySlice.Start; i < ternarySlice.End; i++ {
		s.Set(i, 4)
	}
	rel := crypto.NewLinearRelation(A, s)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel).
		WithTernarySlice(ternarySlice)
	params := fastens20.GeneratePublicParameters(rel, config)
	ExecuteAndTestCorrectness(t, s, params)
}

func TestManyRows(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := bfvRing.D
	n := bfvRing.D
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(n, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(rel, config)
	ExecuteAndTestCorrectness(t, s, params)
}

func TestMultiSplitFull(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := bfvRing.D
	n := bfvRing.D * 4
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(n, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(rel, config)
	ExecuteAndTestCorrectness(t, s, params)
}

func TestAllLevels(t *testing.T) {
	bfvRing := fastmath.BFVFullRing()
	m := 16
	n := bfvRing.D
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(n, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(rel, config)
	ExecuteAndTestCorrectness(t, s, params)
}

func TestAllLevelsShortRing(t *testing.T) {
	// Create the linear relation on the BFV ring.
	bfvRing := fastmath.BFVFullRing()
	m := bfvRing.D
	n := bfvRing.D
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(n, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	// Run the protocol over a smaller ring of degree 2^7.
	commitmentRing := fastmath.BFVFullShortCommtRing(7)
	rebasedRel := rel.Rebase(commitmentRing)
	config := fastens20.DefaultProtocolConfig(commitmentRing, rebasedRel)
	if config.NumSplits() != bfvRing.D/commitmentRing.D {
		t.Errorf("num splits not correct")
	}
	params := fastens20.GeneratePublicParameters(rebasedRel, config)
	ExecuteAndTestCorrectness(t, rebasedRel.S, params)
}

func TestPerformanceSimple(tst *testing.T) {
	bfvRing := fastmath.BFVFullRing()
	m := bfvRing.D
	n := bfvRing.D
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(n, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(rel, config)
	prover := fastens20.NewProver(params)
	verifier := fastens20.NewVerifier(params)
	var t0, t *fastmath.PolyNTTVec
	var w *fastmath.PolyNTTMatrix
	var ps fastens20.ProverState
	{
		start_t := time.Now()
		t0, t, w, ps = prover.CommitToMessage(s)
		end_t := time.Now()
		delta_t := end_t.Sub(start_t)
		tst.Logf("Prover.CommitToMessage execution took %dms", delta_t.Milliseconds())
	}
	var alpha *fastmath.PolyNTTVec
	var gamma *fastmath.IntMatrix
	var vs fastens20.VerifierState
	{
		start_t := time.Now()
		alpha, gamma, vs = verifier.CreateMasks(t0, t, w)
		end_t := time.Now()
		delta_t := end_t.Sub(start_t)
		tst.Logf("Verifier.CreateMasks execution took %dms", delta_t.Milliseconds())
	}
	var h, v *fastmath.PolyNTT
	var vp *fastmath.PolyNTTVec
	{
		start_t := time.Now()
		t, h, v, vp, ps = prover.CommitToRelation(alpha, gamma, ps)
		end_t := time.Now()
		delta_t := end_t.Sub(start_t)
		tst.Logf("Prover.CommitToRelation execution took %dms", delta_t.Milliseconds())
	}
	var c *fastmath.Poly
	{
		start_t := time.Now()
		c, vs = verifier.CreateChallenge(t, h, v, vp, vs)
		end_t := time.Now()
		delta_t := end_t.Sub(start_t)
		tst.Logf("Verifier.CreateChallenge execution took %dms", delta_t.Milliseconds())
	}
	var z *fastmath.PolyNTTMatrix
	{
		var err error
		start_t := time.Now()
		z, ps, err = prover.MaskedOpening(c, ps)
		for err != nil {
			z, ps, err = prover.MaskedOpening(c, ps)
		}
		end_t := time.Now()
		delta_t := end_t.Sub(start_t)
		tst.Logf("Prover.MaskedOpening execution took %dms", delta_t.Milliseconds())
	}
	{
		start_t := time.Now()
		if !verifier.Verify(z, vs) {
			tst.Logf("verification failed")
		}
		end_t := time.Now()
		delta_t := end_t.Sub(start_t)
		tst.Logf("Verifier.Verify execution took %dms", delta_t.Milliseconds())
	}
}
