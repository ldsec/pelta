package ens20

import "github.com/ldsec/codeBase/commitment/math"

type VerifierState struct {
	sig   math.Automorphism
	t0    math.Vector
	t     math.Vector
	w     math.Matrix
	alpha math.Vector
	gamma math.Matrix
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
	alpha := NewRandomPolynomialVector(vf.settings.k*vf.settings.numSplits, vf.settings.baseRing, vf.settings.uniformSampler)
	gamma := NewRandomIntegerMatrix(vf.settings.k, vf.settings.m, vf.settings.q)
	// Create the automorphism.
	sig := math.NewAutomorphism(int64(vf.settings.d), int64(vf.settings.k))
	return alpha, gamma, VerifierState{t0: t0, t: t, w: w, alpha: alpha, gamma: gamma, sig: sig}

}

// CreateChallenge returns a challenge for the relation commitment.
// Returns c
func (vf Verifier) CreateChallenge(t math.Vector, h, v math.Polynomial, vp math.Vector, state VerifierState) (math.Polynomial, VerifierState) {
	// Verifier challenge generation.
	c := NewRandomPolynomial(vf.settings.baseRing, vf.settings.uniformSampler)
	// Update the state.
	state.t = t
	state.c = c
	state.h = h
	state.v = v
	state.vp = vp
	return c, state
}

// Verify verifies a masked opening for the commitment to the message, verifying the proof.
// Returns true iff the proof is valid
func (vf Verifier) Verify(z math.Matrix, state VerifierState) bool {
	maskedOpeningTest := z.AllRows(
		func(zi math.Vector, i int) bool {
			rhs := state.w.Row(i).Copy().Add(state.t0.Mul(state.sig.Permute(int64(i), state.c)))
			return zi.AsCoeffs().L2Norm() < vf.settings.beta && vf.publicParams.B0.MulVec(zi).Eq(rhs)
		})
	if !maskedOpeningTest {
		return false
	}
	// Constructing f
	f := math.NewMatrixFromDimensions(vf.settings.k, vf.settings.numSplits).Populate(
		func(i int, j int) math.RingElement {
			return vf.publicParams.b.Row(j).Dot(z.Row(i)).Add(state.t.Element(j).Copy().Mul(state.sig.Permute(int64(i), state.c)).Neg())
		})
	f2 := vf.publicParams.b.Row(vf.settings.numSplits + 2).Dot(z.Row(0)).Add(state.c.Mul(state.t.Element(vf.settings.numSplits + 2)).Neg())
	f3 := vf.publicParams.b.Row(vf.settings.numSplits + 3).Dot(z.Row(0)).Add(state.c.Mul(state.t.Element(vf.settings.numSplits + 3)).Neg())
	vTest := math.NewMatrixFromDimensions(vf.settings.k, vf.settings.numSplits).Populate(
		func(i int, j int) math.RingElement {
			p1 := f.Element(i, j).Copy()
			p2 := f.Element(i, j).Copy().Add(state.sig.Permute(int64(i), state.c))
			p3 := f.Element(i, j).Copy().Add(state.sig.Permute(int64(i), state.c).Neg())
			return state.alpha.Element(vf.settings.numSplits*i + j).Mul(state.sig.Permute(int64(-i), p1.Mul(p2).Mul(p3).(math.Polynomial)))
		}).Sum().
		Add(f2).
		Add(state.c.Mul(f3)).
		Eq(state.v)
	if !vTest {
		return false
	}
	hTest := math.NewVectorFromSize(vf.settings.k).All(
		func(_ math.RingElement, i int) bool {
			return state.h.Coeff(i) == 0
		})
	if !hTest {
		return false
	}
	// Reconstruct psi
	At := vf.publicParams.A.Copy().AsMatrix().Transpose()
	psi := math.NewMatrixFromDimensions(vf.settings.k, vf.settings.numSplits).PopulateRows(
		func(mu int) math.Vector {
			tmp := At.Copy().AsMatrix().MulVec(state.gamma.Row(mu))
			return SplitInvNTT(tmp, vf.settings.numSplits, vf.settings.d, vf.settings.baseRing)
		})
	// Reconstruct the commitment to f
	invk := math.NewModInt(int64(vf.settings.k), vf.settings.q).Inv()
	tao := math.NewVectorFromSize(vf.settings.k).Populate(
		func(mu int) math.RingElement {
			tmp := math.NewVectorFromSize(vf.settings.k).Populate(
				func(v int) math.RingElement {
					tmp2 := math.NewVectorFromSize(vf.settings.numSplits).Populate(
						func(j int) math.RingElement {
							dec := math.NewOnePolynomial(vf.settings.baseRing).
								Scale(vf.publicParams.u.Dot(state.gamma.Row(mu)).(*math.ModInt).Uint64()).
								Neg()
							return psi.Element(mu, j).Copy().
								Mul(state.t.Element(j)).
								Scale(uint64(vf.settings.d)).
								Add(dec)
						}).Sum()
					return state.sig.Permute(int64(v), tmp2.(math.Polynomial))
				}).Sum().(math.Polynomial)
			return Lmu(mu, tmp, invk)
		}).Sum()
	// Verify the commitments
	testResult := math.NewVectorFromSize(vf.settings.k).Populate(
		func(i int) math.RingElement {
			return math.NewVectorFromSize(vf.settings.k).Populate(
				func(mu int) math.RingElement {
					tmp := math.NewVectorFromSize(vf.settings.k).Populate(
						func(v int) math.RingElement {
							tmp2 := math.NewVectorFromSize(vf.settings.numSplits).Populate(
								func(j int) math.RingElement {
									index := (i - v) % vf.settings.k
									return psi.Element(mu, j).Copy().
										Scale(uint64(vf.settings.d)).
										Mul(vf.publicParams.b.Row(j).Copy().AsVec().Dot(z.Row(index)))
								}).Sum()
							return state.sig.Permute(int64(v), tmp2.(math.Polynomial))
						}).Sum().(math.Polynomial)
					return Lmu(mu, tmp, invk)
				}).Sum().Add(vf.publicParams.b.Row(vf.settings.numSplits + 1).Dot(z.Row(i)))
		}).All(
		func(lhs math.RingElement, i int) bool {
			rhsAdd := state.sig.Permute(int64(i), state.c).Mul(
				tao.Copy().Add(state.t.Element(vf.settings.numSplits + 1)).Copy().Add(state.h.Neg()))
			rhs := state.vp.Element(i).Copy().Add(rhsAdd)
			return lhs.Eq(rhs)
		})
	return testResult
}
