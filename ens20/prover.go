package ens20

import (
	"github.com/ldsec/codeBase/commitment/math"
)

type ProverState struct {
	g    math.Polynomial
	sHat math.Vector
	sig  math.Automorphism
	r    math.Vector
	t0   math.Vector
	t    math.Vector
	y    math.Matrix
	w    math.Matrix
	v    math.Polynomial
	psi  math.Matrix
	h    math.Polynomial
	vp   math.Vector
	z    math.Matrix
}

type Prover struct {
	settings     Settings
	publicParams PublicParams
}

func NewDummyPublicParameters(s math.Vector, settings Settings) PublicParams {
	A := NewRandomIntegerMatrix(settings.m, settings.n, settings.q)
	// As = u
	u := A.Copy().AsMatrix().MulVec(s)
	bSize := settings.numSplits + 3
	B0 := NewRandomPolynomialMatrix(settings.kappa, settings.lambda+settings.kappa+bSize, settings.baseRing, settings.uniformSampler)
	b := NewRandomPolynomialMatrix(bSize, B0.Cols(), settings.baseRing, settings.uniformSampler)

	return PublicParams{
		A:  A,
		u:  u,
		B0: B0,
		b:  b,
	}
}

func NewProver(publicParams PublicParams, settings Settings) Prover {
	return Prover{settings, publicParams}
}

// CommitToMessage commits to the given secret message s.
// Returns t0, t, w
func (p Prover) CommitToMessage(s math.Vector) (math.Vector, math.Vector, math.Matrix, ProverState) {
	// Split the message into polynomial space.
	sHat := SplitInvNTT(s, p.settings.numSplits, p.settings.d, p.settings.baseRing)
	bSize := p.settings.numSplits + 3
	// Sample a polynomial g s.t. g_0=...=g_(k-1)=0
	g := math.NewZeroPolynomial(p.settings.baseRing)
	p.settings.uniformSampler.Read(g.Ref)
	for i := 0; i < p.settings.k; i++ {
		g.SetCoefficient(i, 0)
	}
	// Sample the randomness.
	r := NewRandomPolynomialVector(p.publicParams.B0.Cols()*p.settings.d, p.settings.baseRing, p.settings.ternarySampler)
	// Compute the commitments.
	t0 := p.publicParams.B0.Copy().AsMatrix().MulVec(r)
	t := math.NewVectorFromSize(bSize).Populate(func(i int) math.RingElement {
		if i > p.settings.numSplits {
			return math.NewZeroPolynomial(p.settings.baseRing)
		}
		res := p.publicParams.b.Row(i).Dot(r)
		if i == p.settings.numSplits {
			return res.Add(g)
		} else {
			return res.Add(sHat.Element(i))
		}
	})
	// Create the masks.
	y := NewRandomPolynomialMatrix(p.settings.k, r.Length(), p.settings.baseRing, p.settings.gaussianSampler)
	w := math.NewMatrixFromDimensions(p.settings.k, r.Length()).PopulateRows(func(i int) math.Vector {
		return p.publicParams.B0.Copy().AsMatrix().MulVec(y.Row(i))
	})
	// Create the automorphism.
	sig := math.NewAutomorphism(int64(p.settings.d), int64(p.settings.k))
	return t0, t, w, ProverState{sHat: sHat, g: g, r: r, t0: t0, sig: sig, t: t, y: y, w: w}
}

// CommitToRelation commits to the ternary structure of the secret and the knowledge of the secret s, s.t. As = u.
// Returns t, h, v, vp
func (p Prover) CommitToRelation(alpha math.Vector, gamma math.Matrix, state ProverState) (math.Vector, math.Polynomial, math.Polynomial, math.Vector, ProverState) {
	// Prover further set up.
	sum1 := math.NewMatrixFromDimensions(p.settings.k, p.settings.numSplits).
		Populate(func(i int, j int) math.RingElement {
			tmp := state.sig.Permute(-1,
				state.sHat.Element(j).Copy().(math.Polynomial).
					Scale(uint64(3)).
					Mul(p.publicParams.b.Row(j).Dot(state.y.Row(i)).Pow(2)).(math.Polynomial))
			return alpha.Element(i*p.settings.numSplits + j).Copy().Mul(tmp)
		}).Sum().Neg()
	sum2 := math.NewMatrixFromDimensions(p.settings.k, p.settings.numSplits).
		Populate(func(i int, j int) math.RingElement {
			tmp := state.sig.Permute(-1,
				state.sHat.Element(j).Copy().(math.Polynomial).
					Pow(2).(math.Polynomial).
					Scale(3).
					Add(math.NewZeroPolynomial(p.settings.baseRing).One().Neg()).
					Mul(p.publicParams.b.Row(j).Dot(state.y.Row(i)).Pow(2)).(math.Polynomial))
			return alpha.Element(i*p.settings.numSplits + j).Copy().Mul(tmp)
		}).Sum()
	sum3 := math.NewMatrixFromDimensions(p.settings.k, p.settings.numSplits).
		Populate(func(i int, j int) math.RingElement {
			tmp := state.sig.Permute(-1,
				p.publicParams.b.Row(j).
					Dot(state.y.Row(i)).
					Pow(3).(math.Polynomial))
			return alpha.Element(i*p.settings.numSplits + j).Copy().Mul(tmp)
		}).Sum()
	state.t.SetElementAtIndex(p.settings.numSplits+1,
		p.publicParams.b.Row(p.settings.numSplits+1).Copy().AsVec().
			Dot(state.r).
			Add(p.publicParams.b.Row(p.settings.numSplits+2).Copy().AsVec().Dot(state.y.Row(0))).
			Add(sum1))
	state.t.SetElementAtIndex(p.settings.numSplits+2,
		p.publicParams.b.Row(p.settings.numSplits+2).Copy().AsVec().
			Dot(state.r).
			Add(sum2))
	v := p.publicParams.b.Row(p.settings.numSplits + 1).Copy().AsVec().Dot(state.y.Row(0)).Add(sum3).(math.Polynomial)
	At := p.publicParams.A.Copy().AsMatrix().Transpose()
	psi := math.NewMatrixFromDimensions(p.settings.k, p.settings.numSplits).PopulateRows(
		func(mu int) math.Vector {
			tmp := At.Copy().AsMatrix().MulVec(gamma.Row(mu))
			return SplitInvNTT(tmp, p.settings.numSplits, p.settings.d, p.settings.baseRing)
		})
	invk := math.NewModInt(int64(p.settings.k), p.settings.q).Inv()
	gMask := math.NewVectorFromSize(p.settings.k).Populate(
		func(mu int) math.RingElement {
			tmp := math.NewVectorFromSize(p.settings.k).Populate(
				func(v int) math.RingElement {
					dec := math.NewZeroPolynomial(p.settings.baseRing).One().(math.Polynomial).Scale(
						p.publicParams.u.Copy().AsVec().Dot(gamma.Row(mu)).(*math.ModInt).Uint64())
					return state.sig.Permute(int64(v),
						math.NewVectorFromSize(p.settings.numSplits).Populate(
							func(j int) math.RingElement {
								return psi.Element(mu, j).Copy().
									Mul(state.sHat.Element(j)).(math.Polynomial).
									Scale(uint64(p.settings.d))
							}).
							Sum().
							Add(dec.Neg()).(math.Polynomial))
				}).Sum().(math.Polynomial)
			return tmp.RRot(mu).Mul(invk)
		}).Sum()
	h := state.g.Copy().Add(gMask).(math.Polynomial)
	vp := math.NewVectorFromSize(p.settings.k).Populate(
		func(i int) math.RingElement {
			return math.NewVectorFromSize(p.settings.k).Populate(
				func(mu int) math.RingElement {
					tmp2 := math.NewMatrixFromDimensions(p.settings.k, p.settings.numSplits).Populate(
						func(v, j int) math.RingElement {
							return state.sig.Permute(int64(v),
								p.publicParams.b.Row(j).
									Mul(psi.Element(mu, j)).
									Scale(uint64(p.settings.d)).AsVec().
									Dot(state.y.Row(i-v)).(math.Polynomial))
						}).Sum().(math.Polynomial)
					return tmp2.RRot(mu).Mul(invk)
				}).
				Sum().
				Add(p.publicParams.b.Row(p.settings.numSplits).Dot(state.y.Row(i)))
		})
	// Update the state.
	state.v = v
	state.psi = psi
	state.h = h
	state.vp = vp
	return state.t, h, v, vp, state
}

// MaskedOpening returns the masked openings to the commitments.
// Returns z
func (p Prover) MaskedOpening(c math.Polynomial, state ProverState) (math.Matrix, ProverState, error) {
	// Masked openings.
	z := math.NewMatrixFromDimensions(p.settings.k, state.r.Length()).PopulateRows(
		func(i int) math.Vector {
			sigc := state.sig.Permute(int64(i), c.Copy().(math.Polynomial))
			return state.y.Row(i).Copy().Add(state.r.Mul(sigc)).AsVec()
		})
	// TODO: check the infinity norm
	_ = z.AllRows(
		func(zi math.Vector, i int) bool {
			return true
		})
	// Update the state.
	state.z = z
	return z, state, nil
}
