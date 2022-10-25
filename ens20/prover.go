package ens20

import (
	"github.com/ldsec/codeBase/commitment/math"
)

type ProverState struct {
	G    math.Polynomial
	SHat math.PolyVector
	Sig  math.Automorphism
	R    math.Vector
	T0   math.Vector
	T    math.Vector
	Y    math.Matrix
	W    math.Matrix
	V    math.Polynomial
	Psi  math.Matrix
	H    math.Polynomial
	Vp   math.Vector
	Z    math.Matrix
}

type Prover struct {
	settings     Settings
	publicParams PublicParams
}

func NewProver(publicParams PublicParams, settings Settings) Prover {
	return Prover{settings, publicParams}
}

// CommitToMessage commits to the given secret message s.
// Returns t0, t, w
func (p Prover) CommitToMessage(s math.IntVector) (math.Vector, math.Vector, math.Matrix, ProverState) {
	// Split the message into polynomial space.
	sHat := SplitInvNTT(s, p.settings.NumSplits, p.settings.BaseRing)
	// Sample a polynomial g s.t. g_0=...=g_{K-1}=0
	g := NewRandomPolynomial(p.settings.BaseRing, p.settings.UniformSampler)
	for i := 0; i < p.settings.K; i++ {
		g.SetCoefficient(i, 0)
	}
	// Sample the randomness.
	r := NewRandomPolynomialVector(p.publicParams.B0.Cols(), p.settings.BaseRing, p.settings.TernarySampler)
	// Compute the commitments.
	t0 := p.publicParams.B0.Copy().AsMatrix().MulVec(r)
	t := math.NewVectorFromSize(p.settings.NumSplits + 3).Populate(
		func(j int) math.RingElement {
			if j >= p.settings.NumSplits {
				return math.NewZeroPolynomial(p.settings.BaseRing)
			}
			return p.publicParams.B.Row(j).Copy().AsVec().Dot(r).Add(sHat.Element(j))
		})
	// t[n/d] = b[n/d] * r + g
	t.SetElementAtIndex(p.settings.NumSplits, p.publicParams.B.Row(p.settings.NumSplits).Copy().AsVec().Dot(r).Add(g))
	// Create the masks.
	y := NewRandomPolynomialMatrix(p.settings.K, r.Length(), p.settings.BaseRing, p.settings.GaussianSampler)
	w := math.NewMatrixFromDimensions(p.settings.K, p.settings.Kappa).PopulateRows(
		func(i int) math.Vector {
			return p.publicParams.B0.Copy().AsMatrix().MulVec(y.Row(i))
		})
	// Create the automorphism.
	sig := math.NewAutomorphism(int64(p.settings.D), int64(p.settings.K))
	return t0, t, w, ProverState{SHat: sHat, G: g, R: r, T0: t0, Sig: sig, T: t, Y: y, W: w}
}

// CommitToRelation commits to the ternary structure of the secret and the knowledge of the secret s, s.t. As = U.
// Returns t, h, v, vp
func (p Prover) CommitToRelation(alpha math.Vector, gamma math.Matrix, state ProverState) (math.Vector, math.Polynomial, math.Polynomial, math.Vector, ProverState) {
	// Prover further set up.
	sum1 := CommitmentSum(p.settings.K, p.settings.NumSplits, alpha.AsPolyVec(), state.Sig,
		func(i int, j int) math.Polynomial {
			// (b[j] * y[i])^2
			tmp := p.publicParams.B.Row(j).Copy().AsVec().Dot(state.Y.Row(i)).Pow(2)
			// 3s[j] * (b[j] * y[i])^2
			return state.SHat.Element(j).Copy().(math.Polynomial).
				Scale(uint64(3)).
				Mul(tmp).(math.Polynomial)
		})
	sum2 := CommitmentSum(p.settings.K, p.settings.NumSplits, alpha.AsPolyVec(), state.Sig,
		func(i int, j int) math.Polynomial {
			// b[j] * y[i]
			tmp := p.publicParams.B.Row(j).Copy().AsVec().Dot(state.Y.Row(i))
			// 3s[j]^2
			tmp2 := state.SHat.Element(j).Copy().(math.Polynomial).
				Pow(2).(math.Polynomial).
				Scale(3)
			// (3s[j]^2-1) (b[j] * y[i])
			return tmp2.
				Sub(math.NewOnePolynomial(p.settings.BaseRing)).
				Mul(tmp).(math.Polynomial)
		})
	sum3 := CommitmentSum(p.settings.K, p.settings.NumSplits, alpha.AsPolyVec(), state.Sig,
		func(i int, j int) math.Polynomial {
			// (b[j] * y[i])^3
			return p.publicParams.B.Row(j).Copy().AsVec().
				Dot(state.Y.Row(i)).
				Pow(3).(math.Polynomial)
		})
	state.T.SetElementAtIndex(p.settings.NumSplits+1,
		p.publicParams.B.Row(p.settings.NumSplits+1).Copy().AsVec().
			Dot(state.R).
			Add(p.publicParams.B.Row(p.settings.NumSplits+2).Copy().AsVec().Dot(state.Y.Row(0))).
			Sub(sum1))
	state.T.SetElementAtIndex(p.settings.NumSplits+2,
		p.publicParams.B.Row(p.settings.NumSplits+2).Copy().AsVec().
			Dot(state.R).
			Add(sum2))
	v := p.publicParams.B.Row(p.settings.NumSplits + 1).Copy().AsVec().
		Dot(state.Y.Row(0)).
		Add(sum3).(math.Polynomial)
	At := p.publicParams.A.Copy().AsMatrix().Transpose()
	psi := math.NewMatrixFromDimensions(p.settings.K, p.settings.NumSplits).PopulateRows(
		func(mu int) math.Vector {
			tmp := At.Copy().AsMatrix().MulVec(gamma.Row(mu)).AsIntVec()
			return SplitInvNTT(tmp, p.settings.NumSplits, p.settings.BaseRing).AsVec()
		})
	invk := math.NewModInt(int64(p.settings.K), p.settings.Q).Inv().Uint64()
	gMask := LmuSum(p.settings.K, invk, state.Sig,
		func(mu int, v int) math.Polynomial {
			// (u * gamma_mu)
			mul := p.publicParams.U.Copy().AsVec().Dot(gamma.Row(mu)).(*math.ModInt).Uint64()
			dec := math.NewOnePolynomial(p.settings.BaseRing).Scale(mul)
			sum := math.NewVectorFromSize(p.settings.NumSplits).Populate(
				func(j int) math.RingElement {
					return psi.Element(mu, j).Copy().
						Mul(state.SHat.Element(j)).(math.Polynomial).
						Scale(uint64(p.settings.D))
				}).Sum()
			return sum.Sub(dec).(math.Polynomial)
		})
	h := state.G.Copy().Add(gMask).(math.Polynomial)
	vp := math.NewVectorFromSize(p.settings.K).Populate(
		func(i int) math.RingElement {
			// b[n/d] * y[i]
			add := p.publicParams.B.Row(p.settings.NumSplits).Copy().AsVec().Dot(state.Y.Row(i))
			outerSum := LmuSumOuter(p.settings.K, p.settings.NumSplits, invk, state.Sig,
				func(mu int, v int, j int) math.Polynomial {
					index := math.Mod(i-v, p.settings.K)
					return p.publicParams.B.Row(j).Copy().
						Mul(psi.Element(mu, j)).
						Scale(uint64(p.settings.D)).AsVec().
						Dot(state.Y.Row(index)).(math.Polynomial)
				})
			return outerSum.Add(add)
		})
	// Update the state.
	state.V = v
	state.Psi = psi
	state.H = h
	state.Vp = vp
	return state.T, h, v, vp, state
}

// MaskedOpening returns the masked openings to the commitments.
// Returns z
func (p Prover) MaskedOpening(c math.Polynomial, state ProverState) (math.Matrix, ProverState, error) {
	// Masked openings.
	z := math.NewMatrixFromDimensions(p.settings.K, state.R.Length()).PopulateRows(
		func(i int) math.Vector {
			sigc := state.Sig.Permute(int64(i), c)
			return state.Y.Row(i).Copy().Add(state.R.Copy().Mul(sigc)).AsVec()
		})

	//normInBounds := z.AllRows(
	//	func(zi math.Vector, i int) bool {
	//		infNorm := zi.AsPolyVec().InfNorm(p.settings.Q)
	//		fmt.Println(infNorm)
	//		return float64(infNorm) < p.settings.Beta
	//	})
	//if !normInBounds {
	//	return math.Matrix{}, state, errors.New("infinity norm exceeds the requested bound")
	//}
	// Update the state.
	state.Z = z
	return z, state, nil
}