package tests

import (
	"fmt"
	"math/big"
	"testing"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastens20"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func createRandomRelation(m, n int, relRing fastmath.RingParams) crypto.LinearRelation {
	uni, ter, _ := crypto.GetSamplers(relRing, 128)
	A := fastmath.NewRandomIntMatrixFast(m, n, uni, relRing.BaseRing)
	s := fastmath.NewRandomIntVecFast(n, ter, relRing.BaseRing)
	return crypto.NewLinearRelation(A, s)
}

func executeAndTestCorrectness(tst *testing.T, s *fastmath.IntVec, params fastens20.PublicParams) {
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
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 16
	n := bfvRing.D
	rel := createRandomRelation(m, n, bfvRing)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(config)
	ACopy := params.A.Copy()
	BCopy := params.B.Copy()
	B0Copy := params.B0.Copy()
	prover := fastens20.NewProver(params)
	verifier := fastens20.NewVerifier(params)
	// Commit to the message.
	t0, t, w, ps := prover.CommitToMessage(rel.S)
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
				Pow(2)
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
	rel := createRandomRelation(m, n, bfvRing)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(config)
	executeAndTestCorrectness(t, rel.S, params)
}

func TestMultiReplication(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 16
	n := bfvRing.D
	rel := createRandomRelation(m, n, bfvRing)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel).
		WithReplication(4)
	params := fastens20.GeneratePublicParameters(config)
	executeAndTestCorrectness(t, rel.S, params)
}

func TestMultiSplit(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 16
	n := bfvRing.D * 4
	rel := createRandomRelation(m, n, bfvRing)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(config)
	executeAndTestCorrectness(t, rel.S, params)
}

func TestMultiSplitMultiReplication(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 16
	n := bfvRing.D * 4
	rel := createRandomRelation(m, n, bfvRing)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel).
		WithReplication(4)
	params := fastens20.GeneratePublicParameters(config)
	executeAndTestCorrectness(t, rel.S, params)
}

func TestTernarySlice(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := 16
	n := bfvRing.D * 4
	A := fastmath.NewRandomIntMatrix(m, n, bfvRing.Q, bfvRing.BaseRing)
	s := fastmath.NewRandomIntVec(n, big.NewInt(5), bfvRing.BaseRing)
	ternarySlice := fastmath.NewSlice(n/2, n)
	for i := ternarySlice.Start; i < ternarySlice.End; i++ {
		s.SetForce(i, 4)
	}
	rel := crypto.NewLinearRelation(A, s)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel).
		WithTernarySlice(ternarySlice)
	params := fastens20.GeneratePublicParameters(config)
	executeAndTestCorrectness(t, s, params)
}

func TestDRows(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := bfvRing.D
	n := bfvRing.D
	rel := createRandomRelation(m, n, bfvRing)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(config)
	executeAndTestCorrectness(t, rel.S, params)
}

func TestMultiSplitDRows(t *testing.T) {
	bfvRing := fastmath.BFVZeroLevelRing()
	m := bfvRing.D
	n := bfvRing.D * 4
	rel := createRandomRelation(m, n, bfvRing)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(config)
	executeAndTestCorrectness(t, rel.S, params)
}

func TestFullRing(t *testing.T) {
	bfvRing := fastmath.BFVFullRing()
	m := 16
	n := bfvRing.D
	rel := createRandomRelation(m, n, bfvRing)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(config)
	executeAndTestCorrectness(t, rel.S, params)
}

func TestFullRingDRows(t *testing.T) {
	bfvRing := fastmath.BFVFullRing()
	m := bfvRing.D
	n := bfvRing.D
	rel := createRandomRelation(m, n, bfvRing)
	config := fastens20.DefaultProtocolConfig(bfvRing, rel)
	params := fastens20.GeneratePublicParameters(config)
	executeAndTestCorrectness(t, rel.S, params)
}

func TestFullRingRebaseDRows(t *testing.T) {
	// Create the linear relation on the BFV ring.
	bfvRing := fastmath.BFVFullRing()
	m := bfvRing.D
	n := bfvRing.D
	rel := createRandomRelation(m, n, bfvRing)
	// Run the protocol over a smaller ring of degree 2^7.
	commitmentRing := fastmath.BFVFullShortCommtRing(7)
	rebasedRel := rel.Rebased(commitmentRing)
	config := fastens20.DefaultProtocolConfig(commitmentRing, rebasedRel)
	if config.NumSplits() != bfvRing.D/commitmentRing.D {
		t.Errorf("num splits not correct")
	}
	params := fastens20.GeneratePublicParameters(config)
	executeAndTestCorrectness(t, rebasedRel.S, params)
}

func TestPerformanceBig(tst *testing.T) {
	config := crypto.GetDefaultCryptoConfig()
	bfvRing := config.RingParams
	m := 2 * bfvRing.D
	n := 2 * bfvRing.D
	matrix_name := fmt.Sprintf("performance_test_A_%d_%d.test", m, n)
	vec_name := fmt.Sprintf("performance_test_s_%d.test", n)

	// get the samplers
	uni, ter, _ := crypto.GetSamplers(bfvRing, 128)

	// create the inputs
	A := fastmath.PersistentIntMatrix(matrix_name, func() *fastmath.IntMatrix {
		return fastmath.NewRandomIntMatrixFast(m, n, uni, bfvRing.BaseRing)
	}, bfvRing.BaseRing)
	s := fastmath.PersistentIntVec(vec_name, func() *fastmath.IntVec {
		return fastmath.NewRandomIntVecFast(n, ter, bfvRing.BaseRing)
	}, bfvRing.BaseRing)
	rel := crypto.NewLinearRelation(A, s)

	commitmentRing := fastmath.BFVFullShortCommtRing(7)
	rebasedRel := rel.Rebased(commitmentRing)
	protocolConfig := fastens20.DefaultProtocolConfig(commitmentRing, rebasedRel).
		WithReplication(4)
	params := fastens20.GeneratePublicParameters(protocolConfig)
	if !fastens20.Execute(rebasedRel.S, params) {
		tst.Errorf("execution failed")
	}
}
