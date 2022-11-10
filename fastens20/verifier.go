package fastens20

import (
	"fmt"
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

type VerifierState struct {
	Sig   fastmath.Automorphism
	T0    fastmath.PolyVec
	T     fastmath.PolyVec
	W     fastmath.PolyMatrix
	Alpha fastmath.PolyVec
	Gamma fastmath.IntMatrix
	c     fastmath.Poly
	h     fastmath.Poly
	v     fastmath.Poly
	vp    fastmath.PolyVec
}

type Verifier struct {
	settings     Settings
	publicParams PublicParams
}

func NewVerifier(publicParams PublicParams, settings Settings) Verifier {
	return Verifier{settings, publicParams}
}

// CreateMasks returns the relation masks in response to a message commitment.
// Returns alpha, gamma
func (vf Verifier) CreateMasks(t0 fastmath.PolyVec, t fastmath.PolyVec, w fastmath.PolyMatrix) (fastmath.PolyVec, fastmath.IntMatrix, VerifierState) {
	alpha := fastmath.NewRandomPolyVec(vf.settings.K*vf.settings.NumSplits(), vf.settings.UniformSampler, vf.settings.BaseRing)
	gamma := fastmath.NewRandomIntMatrix(vf.settings.K, vf.settings.M, vf.settings.Q, vf.settings.BaseRing)
	// Create the automorphism.
	sig := fastmath.NewAutomorphism(uint64(vf.settings.D), uint64(vf.settings.K))
	return alpha, gamma, VerifierState{T0: t0, T: t, W: w, Alpha: alpha, Gamma: gamma, Sig: sig}
}

// CreateChallenge returns a challenge for the relation commitment.
// Returns c
func (vf Verifier) CreateChallenge(t fastmath.PolyVec, h, v fastmath.Poly, vp fastmath.PolyVec, state VerifierState) (fastmath.Poly, VerifierState) {
	// Verifier challenge generation.
	c := fastmath.NewRandomPoly(vf.settings.TernarySampler, vf.settings.BaseRing)
	// Update the state.
	state.T = t
	state.c = c
	state.h = h
	state.v = v
	state.vp = vp
	return c, state
}

// Verify verifies a masked opening for the commitment to the message, verifying the proof.
// Returns true iff the proof is valid
func (vf Verifier) Verify(z fastmath.PolyMatrix, state VerifierState) bool {
	maskedOpeningTestResult := z.AllRows(
		func(zi fastmath.PolyVec, i int) bool {
			rhs := state.W.Row(i).Copy().Add(state.T0.Copy().MulAll(state.c.Permute(int64(i), state.Sig)))
			//sizeCheck := zi.Copy().AsVec().NewPolyVec().L2Norm(vf.settings.Q) < vf.settings.Beta
			return vf.publicParams.B0.MulVec(&zi).Eq(rhs) //&& sizeCheck
		})
	if !maskedOpeningTestResult {
		fmt.Println("Verifier.Verify: Failed masked opening verification")
		return false
	}
	// Constructing f
	f := fastmath.NewPolyMatrix(vf.settings.K, vf.settings.NumSplits(), vf.settings.BaseRing)
	f.Populate(func(i int, j int) fastmath.Poly {
		// t[j]*sig^i(c)
		tmp := state.T.Get(j).Copy().MulCoeffs(state.c.Permute(int64(i), state.Sig))
		// (b[j] * z[i]) - (t[j]*sig^i(c))
		return *vf.publicParams.B.Row(j).Dot(z.Row(i)).Add(tmp.Neg())
	})
	// (b[n/d+2] * z[0]) - c*t[n/d+2])
	f2 := vf.publicParams.B.Row(vf.settings.NumSplits() + 1).Copy().
		Dot(z.Row(0)).
		Add(state.c.Copy().
			MulCoeffs(state.T.Get(vf.settings.NumSplits() + 1)).
			Neg())
	// (b[n/d+3] * z[0]) - c*t[n/d+3]
	f3 := vf.publicParams.B.Row(vf.settings.NumSplits() + 2).Copy().
		Dot(z.Row(0)).
		Add(state.c.Copy().
			MulCoeffs(state.T.Get(vf.settings.NumSplits() + 2)).
			Neg())
	vTest := CommitmentSum(vf.settings.K, vf.settings.NumSplits(), state.Alpha, state.Sig,
		func(i int, j int) fastmath.Poly {
			// f[i][j]
			p1 := f.Get(i, j).Copy()
			// f[i][j] + sig^i(c)
			p2 := f.Get(i, j).Copy().Add(state.c.Permute(int64(i), state.Sig))
			// f[i][j] + 2sig^i(c)
			p3 := f.Get(i, j).Copy().Add(state.c.Permute(int64(i), state.Sig).Scale(2))
			// f[i][j] * (f[i][j] + sig^i(c)) * (f[i][j] + 2sig^i(c))
			return *p1.MulCoeffs(p2).MulCoeffs(p3)
		}, vf.settings.BaseRing)
	vTest.Add(f2).Add(f3.Copy().MulCoeffs(&state.c))
	vTestResult := vTest.Eq(&state.v)
	if !vTestResult {
		fmt.Println("Verifier.Verify: Failed relation check")
		return false
	}
	hTestResult := fastmath.NewPolyVec(vf.settings.K).All(
		func(_ fastmath.Poly, i int) bool {
			return state.h.Get(i, 0) == 0
		})
	if !hTestResult {
		fmt.Println("Verifier.Verify: Failed zero-coefficient check")
		return false
	}
	// Reconstruct psi
	At := vf.publicParams.A.Transposed()
	psi := fastmath.NewPolyMatrix(vf.settings.K, vf.settings.NumSplits(), vf.settings.BaseRing)
	psi.PopulateRows(func(mu int) fastmath.PolyVec {
		tmp := At.MulVec(state.Gamma.RowView(mu))
		return SplitInvNTT(tmp, vf.settings.BaseRing)
	})
	// Reconstruct the commitment to f
	invk := big.NewInt(0).ModInverse(big.NewInt(int64(vf.settings.K)), vf.settings.Q).Uint64()
	tao := LmuSum(vf.settings.K, invk, state.Sig,
		func(mu int, v int) fastmath.Poly {
			// (u * gamma_mu)
			mul := vf.publicParams.U.Dot(state.Gamma.RowView(mu))
			dec := fastmath.NewOnePoly(mul, vf.settings.BaseRing)
			// \sum_{j=0}^{numSplits-1} (d*psi[mu][j] * (b[j] * z[i - v]))
			sum := fastmath.NewPolyVec(vf.settings.NumSplits(), vf.settings.BaseRing)
			sum.Populate(func(j int) fastmath.Poly {
				return *psi.Get(mu, j).Copy().
					MulCoeffs(state.T.Get(j)).
					Scale(uint64(vf.settings.D))
			})
			return *sum.Sum().Add(dec.Neg())
		}, vf.settings.BaseRing)
	// Verify the function commitment
	functionCommitmentTest := fastmath.NewPolyVec(vf.settings.K, vf.settings.BaseRing)
	functionCommitmentTest.Populate(func(i int) fastmath.Poly {
		// (b[n/d] * z[i])
		add := vf.publicParams.B.Row(vf.settings.NumSplits()).Dot(z.Row(i))
		outerSum := LmuSumOuter(vf.settings.K, vf.settings.NumSplits(), invk, state.Sig,
			func(mu int, v int, j int) fastmath.Poly {
				index := big.NewInt(0).Mod(big.NewInt(int64(i-v)), big.NewInt(int64(vf.settings.K))).Int64()
				// b[j] * z[i - v]
				mul := vf.publicParams.B.Row(j).Copy().Dot(z.Row(int(index)))
				// d * psi[mu][j] * (b[j] * z[i - v])
				return *psi.Get(mu, j).Copy().
					Scale(uint64(vf.settings.D)).
					MulCoeffs(mul)
			}, vf.settings.BaseRing)
		// outerSum + (b[n/d] * z[i])
		return *outerSum.Add(add)
	})
	functionCommitmentTestResult := functionCommitmentTest.All(
		func(lhs fastmath.Element, i int) bool {
			rhsAdd := state.c.Permute(int64(i), state.Sig).MulCoeffs(
				tao.Copy().Add(state.T.Get(vf.settings.NumSplits())).Copy().Add(state.h.Neg()))
			rhs := state.vp.Get(i).Copy().Add(rhsAdd)
			return lhs.Eq(rhs)
		})
	if !functionCommitmentTestResult {
		fmt.Println("Verifier.Verify: Failed function commitment check")
		return false
	}
	return true
}
