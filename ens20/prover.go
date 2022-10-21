package ens20

import (
	"github.com/ldsec/codeBase/commitment/math"
)

type ProverState struct {
	g    math.Polynomial
	sHat math.PolyVector
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

func NewDummyPublicParameters(s math.IntVector, settings Settings) PublicParams {
	A := NewRandomIntegerMatrix(settings.M, settings.N, settings.Q)
	// As = u
	u := A.Copy().AsMatrix().MulVec(s.AsVec()).AsIntVec()
	bSize := settings.NumSplits + 3
	B0 := NewRandomPolynomialMatrix(settings.Kappa, settings.Lambda+settings.Kappa+bSize, settings.BaseRing, settings.UniformSampler)
	b := NewRandomPolynomialMatrix(bSize, B0.Cols(), settings.BaseRing, settings.UniformSampler)

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
func (p Prover) CommitToMessage(s math.IntVector) (math.Vector, math.Vector, math.Matrix, ProverState) {
	// Split the message into polynomial space.
	sHat := SplitInvNTT(s, p.settings.NumSplits, p.settings.D, p.settings.BaseRing)
	bSize := p.settings.NumSplits + 3
	// Sample a polynomial g s.t. g_0=...=g_{K-1}=0
	g := math.NewZeroPolynomial(p.settings.BaseRing)
	p.settings.UniformSampler.Read(g.Ref)
	for i := 0; i < p.settings.K; i++ {
		g.SetCoefficient(i, 0)
	}
	// Sample the randomness.
	r := NewRandomPolynomialVector(p.publicParams.B0.Cols(), p.settings.BaseRing, p.settings.TernarySampler)
	// Compute the commitments.
	t0 := p.publicParams.B0.Copy().AsMatrix().MulVec(r)
	t := math.NewVectorFromSize(bSize).Populate(func(i int) math.RingElement {
		if i > p.settings.NumSplits {
			return math.NewZeroPolynomial(p.settings.BaseRing)
		}
		res := p.publicParams.b.Row(i).Dot(r)
		if i == p.settings.NumSplits {
			return res.Add(g)
		} else {
			return res.Add(sHat.Element(i))
		}
	})
	// Create the masks.
	y := NewRandomPolynomialMatrix(p.settings.K, r.Length(), p.settings.BaseRing, p.settings.GaussianSampler)
	w := math.NewMatrixFromDimensions(p.settings.K, p.settings.Kappa).PopulateRows(
		func(i int) math.Vector {
			v := p.publicParams.B0.Copy().AsMatrix().MulVec(y.Row(i))
			return v
		})
	// Create the automorphism.
	sig := math.NewAutomorphism(int64(p.settings.D), int64(p.settings.K))
	return t0, t, w, ProverState{sHat: sHat, g: g, r: r, t0: t0, sig: sig, t: t, y: y, w: w}
}

// CommitToRelation commits to the ternary structure of the secret and the knowledge of the secret s, s.t. As = u.
// Returns t, h, v, vp
func (p Prover) CommitToRelation(alpha math.Vector, gamma math.Matrix, state ProverState) (math.Vector, math.Polynomial, math.Polynomial, math.Vector, ProverState) {
	// Prover further set up.
	sum1 := math.NewMatrixFromDimensions(p.settings.K, p.settings.NumSplits).
		Populate(func(i int, j int) math.RingElement {
			tmp := state.sig.Permute(-1,
				state.sHat.Element(j).Copy().(math.Polynomial).
					Scale(uint64(3)).
					Mul(p.publicParams.b.Row(j).Dot(state.y.Row(i)).Pow(2)).(math.Polynomial))
			return alpha.Element(i*p.settings.NumSplits + j).Copy().Mul(tmp)
		}).Sum().Neg()
	sum2 := math.NewMatrixFromDimensions(p.settings.K, p.settings.NumSplits).
		Populate(func(i int, j int) math.RingElement {
			tmp := state.sig.Permute(-1,
				state.sHat.Element(j).Copy().(math.Polynomial).
					Pow(2).(math.Polynomial).
					Scale(3).
					Add(math.NewZeroPolynomial(p.settings.BaseRing).One().Neg()).
					Mul(p.publicParams.b.Row(j).Dot(state.y.Row(i)).Pow(2)).(math.Polynomial))
			return alpha.Element(i*p.settings.NumSplits + j).Copy().Mul(tmp)
		}).Sum()
	sum3 := math.NewMatrixFromDimensions(p.settings.K, p.settings.NumSplits).
		Populate(func(i int, j int) math.RingElement {
			tmp := state.sig.Permute(-1,
				p.publicParams.b.Row(j).
					Dot(state.y.Row(i)).
					Pow(3).(math.Polynomial))
			return alpha.Element(i*p.settings.NumSplits + j).Copy().Mul(tmp)
		}).Sum()
	state.t.SetElementAtIndex(p.settings.NumSplits+1,
		p.publicParams.b.Row(p.settings.NumSplits+1).Copy().AsVec().
			Dot(state.r).
			Add(p.publicParams.b.Row(p.settings.NumSplits+2).Copy().AsVec().Dot(state.y.Row(0))).
			Add(sum1))
	state.t.SetElementAtIndex(p.settings.NumSplits+2,
		p.publicParams.b.Row(p.settings.NumSplits+2).Copy().AsVec().
			Dot(state.r).
			Add(sum2))
	v := p.publicParams.b.Row(p.settings.NumSplits + 1).Copy().AsVec().Dot(state.y.Row(0)).Add(sum3).(math.Polynomial)
	At := p.publicParams.A.Copy().AsMatrix().Transpose()
	psi := math.NewMatrixFromDimensions(p.settings.K, p.settings.NumSplits).PopulateRows(
		func(mu int) math.Vector {
			tmp := At.Copy().AsMatrix().MulVec(gamma.Row(mu)).AsIntVec()
			return SplitInvNTT(tmp, p.settings.NumSplits, p.settings.D, p.settings.BaseRing).AsVec()
		})
	invk := math.NewModInt(int64(p.settings.K), p.settings.Q).Inv()
	gMask := math.NewVectorFromSize(p.settings.K).Populate(
		func(mu int) math.RingElement {
			tmp := math.NewVectorFromSize(p.settings.K).Populate(
				func(v int) math.RingElement {
					dec := math.NewZeroPolynomial(p.settings.BaseRing).One().(math.Polynomial).Scale(
						p.publicParams.u.Copy().AsVec().Dot(gamma.Row(mu)).(*math.ModInt).Uint64())
					return state.sig.Permute(int64(v),
						math.NewVectorFromSize(p.settings.NumSplits).Populate(
							func(j int) math.RingElement {
								return psi.Element(mu, j).Copy().
									Mul(state.sHat.Element(j)).(math.Polynomial).
									Scale(uint64(p.settings.D))
							}).
							Sum().
							Add(dec.Neg()).(math.Polynomial))
				}).Sum().(math.Polynomial)
			return Lmu(mu, tmp, invk)
		}).Sum()
	h := state.g.Copy().Add(gMask).(math.Polynomial)
	vp := math.NewVectorFromSize(p.settings.K).Populate(
		func(i int) math.RingElement {
			return math.NewVectorFromSize(p.settings.K).Populate(
				func(mu int) math.RingElement {
					tmp2 := math.NewMatrixFromDimensions(p.settings.K, p.settings.NumSplits).Populate(
						func(v, j int) math.RingElement {
							return state.sig.Permute(int64(v),
								p.publicParams.b.Row(j).
									Mul(psi.Element(mu, j)).
									Scale(uint64(p.settings.D)).AsVec().
									Dot(state.y.Row(i-v)).(math.Polynomial))
						}).Sum().(math.Polynomial)
					return Lmu(mu, tmp2, invk)
				}).
				Sum().
				Add(p.publicParams.b.Row(p.settings.NumSplits).Dot(state.y.Row(i)))
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
	z := math.NewMatrixFromDimensions(p.settings.K, state.r.Length()).PopulateRows(
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
