package tests

import (
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func ExecuteAndTestCorrectness(tst *testing.T, s *fastmath.IntVec, params fastens20.PublicParams) {
	prover := fastens20.NewProver(params)
	verifier := fastens20.NewVerifier(params)
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
		tst.Errorf("verification failed")
	}
}

func TestConsistency(tst *testing.T) {
	config := crypto.GetDefaultConfig()
	config.N *= 2
	// Create a simple SIS problem instance & its solution.
	A := fastmath.NewRandomIntMatrix(config.M, config.N, config.Q, config.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(config.N, config.BaseRing)
	sisProblem := crypto.NewSISProblem(A, s)
	params := fastens20.NewPublicParameters(sisProblem, config)
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
	f.Populate(func(i, j int) fastmath.PolyNTT {
		// t[j]*sig^i(c)
		tmp := t.Get(j).Copy().Mul(c.Permute(int64(i), params.Sig).NTT())
		// (b[j] * z[i]) - (t[j]*sig^i(c))
		return *params.B.Row(j).Dot(z.Row(i)).Add(tmp.Neg())
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
		func(i, j int) fastmath.PolyNTT {
			// (b[j] * y[i])^2
			tmp := params.B.Row(j).Copy().
				Dot(ps.Y.Row(i)).
				Pow(2, config.Q.Uint64())
			// 3s[j] - 3
			tmp2 := ps.SHat.Get(j).Copy().
				Scale(3).
				Add(fastmath.NewOnePoly(3, config.BaseRing).Neg().NTT())
			// (3s[j]-3) (b[j] * y[i])^2
			return *tmp2.Mul(tmp)
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
		func(i, j int) fastmath.PolyNTT {
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
			return *tmp.Mul(tmp1.Mul(tmp2).Add(tmp3))
		}, params)
	check3RHS := params.B.Row(config.NumSplits() + 2).
		Dot(ps.Y.Row(0)).
		Add(check3Sum.Mul(c.Copy().NTT()).Neg())
	if !f3.Eq(check3RHS) {
		tst.Errorf("f[n/d+3] equivalence test failed")
	}

}

func TestSimple(t *testing.T) {
	config := crypto.GetDefaultConfig()
	A := fastmath.NewRandomIntMatrix(config.M, config.N, config.Q, config.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(config.N, config.BaseRing)
	sisProblem := crypto.NewSISProblem(A, s)
	params := fastens20.NewPublicParameters(sisProblem, config)
	ExecuteAndTestCorrectness(t, s, params)
}

func TestMultiReplication(t *testing.T) {
	config := crypto.GetDefaultConfig()
	config.K = 4
	A := fastmath.NewRandomIntMatrix(config.M, config.N, config.Q, config.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(config.N, config.BaseRing)
	sisProblem := crypto.NewSISProblem(A, s)
	params := fastens20.NewPublicParameters(sisProblem, config)
	ExecuteAndTestCorrectness(t, s, params)
}

func TestMultiSplit(t *testing.T) {
	config := crypto.GetDefaultConfig()
	config.N *= 4
	A := fastmath.NewRandomIntMatrix(config.M, config.N, config.Q, config.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(config.N, config.BaseRing)
	sisProblem := crypto.NewSISProblem(A, s)
	params := fastens20.NewPublicParameters(sisProblem, config)
	ExecuteAndTestCorrectness(t, s, params)
}

func TestMultiSplitMultiReplication(t *testing.T) {
	config := crypto.GetDefaultConfig()
	config.N *= 4
	config.K = 4
	A := fastmath.NewRandomIntMatrix(config.M, config.N, config.Q, config.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(config.N, config.BaseRing)
	sisProblem := crypto.NewSISProblem(A, s)
	params := fastens20.NewPublicParameters(sisProblem, config)
	ExecuteAndTestCorrectness(t, s, params)
}

func TestSubTernaryLength(t *testing.T) {
	config := crypto.GetDefaultConfig()
	// Let only half of the secret s have ternary structure.
	config.N *= 4
	config.TernaryLength *= 2
	A := fastmath.NewRandomIntMatrix(config.M, config.N, config.Q, config.BaseRing)
	s := fastmath.NewRandomTernaryIntVec(config.TernaryLength, config.BaseRing)
	s.Append(fastmath.NewRandomIntVec(config.N-config.TernaryLength, config.Q, config.BaseRing))
	sisProblem := crypto.NewSISProblem(A, s)
	params := fastens20.NewPublicParameters(sisProblem, config)
	ExecuteAndTestCorrectness(t, s, params)
}
