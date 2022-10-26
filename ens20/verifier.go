package ens20

import (
	"fmt"
	"github.com/ldsec/codeBase/commitment/math"
)

type VerifierState struct {
	Sig   math.Automorphism
	T0    math.Vector
	T     math.Vector
	W     math.Matrix
	Alpha math.Vector
	Gamma math.Matrix
	c     math.Polynomial
	h     math.Polynomial
	v     math.Polynomial
	vp    math.Vector
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
func (vf Verifier) CreateMasks(t0 math.Vector, t math.Vector, w math.Matrix) (math.Vector, math.Matrix, VerifierState) {
	alpha := NewRandomPolynomialVector(vf.settings.K*vf.settings.NumSplits(), vf.settings.BaseRing, vf.settings.UniformSampler)
	gamma := NewRandomIntegerMatrix(vf.settings.K, vf.settings.M, vf.settings.Q)
	// Create the automorphism.
	sig := math.NewAutomorphism(int64(vf.settings.D), int64(vf.settings.K))
	return alpha, gamma, VerifierState{T0: t0, T: t, W: w, Alpha: alpha, Gamma: gamma, Sig: sig}
}

// CreateChallenge returns a challenge for the relation commitment.
// Returns c
func (vf Verifier) CreateChallenge(t math.Vector, h, v math.Polynomial, vp math.Vector, state VerifierState) (math.Polynomial, VerifierState) {
	// Verifier challenge generation.
	c := NewRandomPolynomial(vf.settings.BaseRing, vf.settings.TernarySampler)
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
func (vf Verifier) Verify(z math.Matrix, state VerifierState) bool {
	maskedOpeningTestResult := z.AllRows(
		func(zi math.Vector, i int) bool {
			rhs := state.W.Row(i).Copy().Add(state.T0.Copy().Mul(state.Sig.Permute(int64(i), state.c)))
			//sizeCheck := zi.Copy().AsVec().AsPolyVec().L2Norm(vf.settings.Q) < vf.settings.Beta
			return vf.publicParams.B0.Copy().AsMatrix().MulVec(zi).Eq(rhs) //&& sizeCheck
		})
	if !maskedOpeningTestResult {
		fmt.Println("Verifier.Verify: Failed masked opening verification")
		return false
	}
	// Constructing f
	f := math.NewMatrixFromDimensions(vf.settings.K, vf.settings.NumSplits()).Populate(
		func(i int, j int) math.RingElement {
			tmp := state.T.Element(j).Copy().Mul(state.Sig.Permute(int64(i), state.c))
			return vf.publicParams.B.Row(j).Copy().AsVec().Dot(z.Row(i)).Sub(tmp)
		})
	f2 := vf.publicParams.B.Row(vf.settings.NumSplits() + 1).Copy().AsVec().Dot(z.Row(0)).Sub(state.c.Copy().Mul(state.T.Element(vf.settings.NumSplits() + 1)))
	f3 := vf.publicParams.B.Row(vf.settings.NumSplits() + 2).Copy().AsVec().Dot(z.Row(0)).Sub(state.c.Copy().Mul(state.T.Element(vf.settings.NumSplits() + 2)))
	vTest := CommitmentSum(vf.settings.K, vf.settings.NumSplits(), state.Alpha.AsPolyVec(), state.Sig,
		func(i int, j int) math.Polynomial {
			p1 := f.Element(i, j).Copy()
			p2 := f.Element(i, j).Copy().Add(state.Sig.Permute(int64(i), state.c))
			p3 := f.Element(i, j).Copy().Sub(state.Sig.Permute(int64(i), state.c))
			return p1.Mul(p2).Mul(p3).(math.Polynomial)
		}).Add(f2).Add(f3.Copy().Mul(state.c))
	vTestResult := vTest.Eq(state.v)
	if !vTestResult {
		fmt.Println("Verifier.Verify: Failed relation check")
		return false
	}
	hTestResult := math.NewVectorFromSize(vf.settings.K).All(
		func(_ math.RingElement, i int) bool {
			return state.h.Coeff(i) == 0
		})
	if !hTestResult {
		fmt.Println("Verifier.Verify: Failed zero-coefficient check")
		return false
	}
	// Reconstruct psi
	At := vf.publicParams.A.Copy().AsMatrix().Transpose()
	psi := math.NewMatrixFromDimensions(vf.settings.K, vf.settings.NumSplits()).PopulateRows(
		func(mu int) math.Vector {
			tmp := At.Copy().AsMatrix().MulVec(state.Gamma.Row(mu)).AsIntVec()
			return SplitInvNTT(tmp, vf.settings.NumSplits(), vf.settings.BaseRing).AsVec()
		})
	// Reconstruct the commitment to f
	invk := math.NewModInt(int64(vf.settings.K), vf.settings.Q).Inv().Uint64()
	tao := LmuSum(vf.settings.K, invk, state.Sig,
		func(mu int, v int) math.Polynomial {
			// (u * gamma_mu)
			mul := vf.publicParams.U.Copy().AsVec().Dot(state.Gamma.Row(mu)).(*math.ModInt).Uint64()
			dec := math.NewOnePolynomial(vf.settings.BaseRing).Scale(mul)
			// \sum_{j=0}^{numSplits-1} (d*psi[mu][j] * (b[j] * z[i - v]))
			sum := math.NewVectorFromSize(vf.settings.NumSplits()).Populate(
				func(j int) math.RingElement {
					return psi.Element(mu, j).Copy().
						Mul(state.T.Element(j)).(math.Polynomial).
						Scale(uint64(vf.settings.D))
				}).Sum()
			return sum.Sub(dec).(math.Polynomial)
		})
	// Verify the function commitment
	functionCommitmentTest := math.NewVectorFromSize(vf.settings.K).Populate(
		func(i int) math.RingElement {
			// (b[n/d] * z[i])
			add := vf.publicParams.B.Row(vf.settings.NumSplits()).Copy().AsVec().Dot(z.Row(i))
			outerSum := LmuSumOuter(vf.settings.K, vf.settings.NumSplits(), invk, state.Sig,
				func(mu int, v int, j int) math.Polynomial {
					index := math.Mod(i-v, vf.settings.K)
					// b[j] * z[i - v]
					mul := vf.publicParams.B.Row(j).Copy().AsVec().Dot(z.Row(index))
					// d * psi[mu][j] * (b[j] * z[i - v])
					return psi.Element(mu, j).Copy().
						Scale(uint64(vf.settings.D)).
						Mul(mul).(math.Polynomial)
				})
			// outerSum + (b[n/d] * z[i])
			return outerSum.Add(add)
		})
	//functionCommitmentTestRHS := math.NewVectorFromSize(vf.settings.K).Populate(
	//	func(i int) math.RingElement {
	//		rhsAdd := state.Sig.Permute(int64(i), state.c).Mul(
	//			tao.Copy().Add(state.T.Element(vf.settings.NumSplits())).Sub(state.h))
	//		return state.vp.Element(i).Copy().Add(rhsAdd)
	//	})
	functionCommitmentTestResult := functionCommitmentTest.All(
		func(lhs math.RingElement, i int) bool {
			rhsAdd := state.Sig.Permute(int64(i), state.c).Mul(
				tao.Copy().Add(state.T.Element(vf.settings.NumSplits())).Sub(state.h))
			rhs := state.vp.Element(i).Copy().Add(rhsAdd)
			return lhs.Eq(rhs)
		})
	if !functionCommitmentTestResult {
		fmt.Println("Verifier.Verify: Failed function commitment check")
		return false
	}
	return true
}
