package fastens20

import (
	"fmt"
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
)

type VerifierState struct {
	T0    *fastmath.PolyNTTVec
	T     *fastmath.PolyNTTVec
	W     *fastmath.PolyNTTMatrix
	Alpha *fastmath.PolyNTTVec
	Gamma *fastmath.IntMatrix
	c     *fastmath.Poly
	h     *fastmath.PolyNTT
	v     *fastmath.PolyNTT
	vp    *fastmath.PolyNTTVec
}

type Verifier struct {
	params PublicParams
}

func NewVerifier(params PublicParams) Verifier {
	return Verifier{params}
}

// CreateMasks returns the relation masks in response to a message commitment.
// Returns alpha, gamma
func (vf Verifier) CreateMasks(t0 *fastmath.PolyNTTVec, t *fastmath.PolyNTTVec, w *fastmath.PolyNTTMatrix) (*fastmath.PolyNTTVec, *fastmath.IntMatrix, VerifierState) {
	e := logging.LogExecStart("Verifier.CreateMasks", "working")
	defer e.LogExecEnd()
	alpha := fastmath.NewRandomPolyVec(vf.params.config.K*vf.params.config.NumSplits(), vf.params.config.UniformSampler, vf.params.config.BaseRing).NTT()
	gamma := fastmath.NewRandomIntMatrix(vf.params.config.K, vf.params.config.M, vf.params.config.Q, vf.params.config.BaseRing)
	return alpha.Copy(), gamma.Copy(), VerifierState{T0: t0, T: t, W: w, Alpha: alpha, Gamma: gamma}
}

// CreateChallenge returns a challenge for the relation commitment.
// Returns c
func (vf Verifier) CreateChallenge(t *fastmath.PolyNTTVec, h, v *fastmath.PolyNTT, vp *fastmath.PolyNTTVec, state VerifierState) (*fastmath.Poly, VerifierState) {
	e := logging.LogExecStart("Verifier.CreateChallenge", "working")
	defer e.LogExecEnd()
	// Verifier challenge generation.
	c := fastmath.NewRandomPoly(vf.params.config.TernarySampler, vf.params.config.BaseRing)
	// Update the state.
	state.T = t
	state.c = c
	state.h = h
	state.v = v
	state.vp = vp
	return c.Copy(), state
}

// Verify verifies a masked opening for the commitment to the message, verifying the proof.
// Returns true iff the proof is valid
func (vf Verifier) Verify(z *fastmath.PolyNTTMatrix, state VerifierState) bool {
	e0 := logging.LogExecStart("Verifier.Verify", "working")
	defer e0.LogExecEnd()
	// Masked opening check
	e := logging.LogExecStart("Verifier.Verify", "masked opening test")
	maskedOpeningTestResult := z.AllRows(
		func(i int, zi *fastmath.PolyNTTVec) bool {
			perm := state.c.Permute(int64(i), vf.params.Sig).NTT()
			rhs := state.T0.Copy().MulAll(perm).Add(state.W.Row(i))
			// TODO masked opening bound check (need a good delta1 parameter)
			//sizeCheck := zi.Copy().AsVec().NewPolyVec().L2Norm(vf.settings.Q) < vf.settings.Beta
			return vf.params.B0.MulVec(zi).Eq(rhs) //&& sizeCheck
		})
	e.LogExecEnd()
	if !maskedOpeningTestResult {
		fmt.Println("verifier failed masked opening verification")
		return false
	}
	// Zero-coefficient check
	e = logging.LogExecStart("Verifier.Verify", "zero coefficient test")
	hTestResult := true
	hInvNTT := state.h.Copy().InvNTT()
	for i := 0; i < vf.params.config.K; i++ {
		if hInvNTT.Get(i, 0) != 0 {
			hTestResult = false
			break
		}
	}
	e.LogExecEnd()
	if !hTestResult {
		fmt.Println(state.h.String())
		fmt.Println("verifier failed zero-coefficient check")
		return false
	}
	// Constructing f
	e = logging.LogExecStart("Verifier.Verify", "relation test (v check)")
	cNTT := state.c.Copy().NTT()
	f := fastmath.NewPolyMatrix(vf.params.config.K, vf.params.config.NumSplits(), vf.params.config.BaseRing).NTT()
	f.Populate(func(i int, j int) *fastmath.PolyNTT {
		// t[j]*sig^i(c)
		tmp := state.c.Permute(int64(i), vf.params.Sig).NTT().Mul(state.T.Get(j))
		// (b[j] * z[i]) - (t[j]*sig^i(c))
		return vf.params.B.Row(j).Dot(z.Row(i)).Add(tmp.Neg())
	})
	// (b[n/d+2] * z[0]) - c*t[n/d+2])
	f2 := vf.params.B.Row(vf.params.config.NumSplits() + 1).
		Dot(z.Row(0)).
		Add(cNTT.Copy().
			Mul(state.T.Get(vf.params.config.NumSplits() + 1)).
			Neg())
	// (b[n/d+3] * z[0]) - c*t[n/d+3]
	f3 := vf.params.B.Row(vf.params.config.NumSplits() + 2).
		Dot(z.Row(0)).
		Add(cNTT.Copy().
			Mul(state.T.Get(vf.params.config.NumSplits() + 2)).
			Neg())
	vTest := CommitmentSum(vf.params.config.K, vf.params.config.NumTernarySplits(), state.Alpha,
		func(i int, j int) *fastmath.PolyNTT {
			// Ternary subvector index
			jp := vf.params.config.TernarySlice.Start/vf.params.config.D + j
			// f[i][j]
			p1 := f.Get(i, jp).Copy()
			// f[i][j] + sig^i(c)
			sigic := state.c.Permute(int64(i), vf.params.Sig).NTT()
			p2 := f.Get(i, jp).Copy().Add(sigic)
			// f[i][j] + 2sig^i(c)
			p3 := p2.Copy().Add(sigic)
			// f[i][j] * (f[i][j] + sig^i(c)) * (f[i][j] + 2sig^i(c))
			return p1.Mul(p2).Mul(p3)
		}, vf.params)
	vTest.Add(f2).Add(f3.Mul(cNTT))
	vTestResult := vTest.Eq(state.v)
	if !vTestResult {
		fmt.Println("verifier failed relation check")
		return false
	}
	e.LogExecEnd()
	e = logging.LogExecStart("Verifier.Verify", "function commitment test")
	// Reconstruct psi
	e2 := logging.LogExecStart("Verifier.Verify", "recalculating psi")
	psi := fastmath.NewPolyMatrix(vf.params.config.K, vf.params.config.NumSplits(), vf.params.config.BaseRing).NTT()
	psi.PopulateRows(func(mu int) *fastmath.PolyNTTVec {
		tmp := vf.params.At.MulVec(state.Gamma.RowView(mu))
		return SplitInvNTT(tmp, vf.params).NTT()
	})
	e2.LogExecEnd()
	// Reconstruct the commitment to f
	e2 = logging.LogExecStart("Verifier.Verify", "calculating tao, the function commitment")
	tao := LmuSum(vf.params.config.K, vf.params.config.InvK,
		func(mu int, v int) *fastmath.PolyNTT {
			// (u * gamma_mu)
			mul := vf.params.U.Dot(state.Gamma.RowView(mu))
			dec := fastmath.NewOnePoly(mul, vf.params.config.BaseRing).NTT()
			// \sum_{j=0}^{numSplits-1} (d*psi[mu][j] * (b[j] * z[i - v]))
			presum := fastmath.NewPolyVec(vf.params.config.NumSplits(), vf.params.config.BaseRing).NTT()
			presum.Populate(func(j int) *fastmath.PolyNTT {
				return psi.Get(mu, j).Copy().
					Mul(state.T.Get(j)).
					Scale(uint64(vf.params.config.D))
			})
			return presum.Sum().Add(dec.Neg())
		}, vf.params)
	e2.LogExecEnd()
	e2 = logging.LogExecStart("Verifier.Verify", "calculating lhs")
	// Verify the function commitment
	functionCommitmentTest := fastmath.NewPolyVec(vf.params.config.K, vf.params.config.BaseRing).NTT()
	functionCommitmentTest.Populate(func(i int) *fastmath.PolyNTT {
		// (b[n/d + 1] * z[i])
		add := vf.params.B.Row(vf.params.config.NumSplits()).Dot(z.Row(i))
		outerSum := LmuSumOuter(vf.params.config.K, vf.params.config.NumSplits(), vf.params.config.InvK,
			func(mu int, v int, j int) *fastmath.PolyNTT {
				index := big.NewInt(0).
					Mod(big.NewInt(int64(i-v)),
						big.NewInt(int64(vf.params.config.K))).
					Int64()
				// b[j] * z[i - v]
				dotResult := vf.params.B.Row(j).Dot(z.Row(int(index)))
				// d * psi[mu][j] * (b[j] * z[i - v])
				return dotResult.Scale(uint64(vf.params.config.D)).
					Mul(psi.Get(mu, j))
			}, vf.params)
		// outerSum + (b[n/d] * z[i])
		return outerSum.Add(add)
	})
	e2.LogExecEnd()
	hNeg := state.h.Copy().Neg()
	functionCommitmentTestResult := functionCommitmentTest.All(
		func(i int, lhs *fastmath.PolyNTT) bool {
			rhsAdd := state.c.Permute(int64(i), vf.params.Sig).NTT().Mul(
				tao.Copy().Add(state.T.Get(vf.params.config.NumSplits())).Add(hNeg))
			rhs := state.vp.Get(i).Add(rhsAdd)
			return lhs.Eq(rhs)
		})
	e.LogExecEnd()
	if !functionCommitmentTestResult {
		fmt.Println("verifier failed function commitment check")
		return false
	}
	return true
}
