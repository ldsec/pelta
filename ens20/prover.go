package ens20

import (
	"github.com/ldsec/codeBase/commitment/math"
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
)

type ProverState struct {
	G    rings.Polynomial
	SHat rings.PolyVector
	Sig  math.Automorphism
	R    algebra.Vector
	T0   algebra.Vector
	T    algebra.Vector
	Y    algebra.Matrix
	W    algebra.Matrix
	V    rings.Polynomial
	Psi  algebra.Matrix
	H    rings.Polynomial
	Vp   algebra.Vector
	Z    algebra.Matrix
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
func (p Prover) CommitToMessage(s rings.ZIntVector) (algebra.Vector, algebra.Vector, algebra.Matrix, ProverState) {
	// Rebase the message into polynomial space.
	sHat := SplitInvNTT(s, p.settings.NumSplits(), p.settings.Q, p.settings.BaseRing)
	// Sample a polynomial g s.t. g_0=...=g_{K-1}=0
	g := math.NewRandomPolynomial(p.settings.BaseRing, p.settings.UniformSampler)
	for i := 0; i < p.settings.K; i++ {
		g.SetCoefficient(i, 0)
	}
	// Sample the randomness.
	r := math.NewRandomPolynomialVector(p.publicParams.B0.Cols(), p.settings.BaseRing, p.settings.TernarySampler)
	// Compute the commitments.
	t0 := p.publicParams.B0.Copy().AsMatrix().MulVec(r)
	t := algebra.NewVectorFromSize(p.settings.NumSplits() + 3).Populate(
		func(j int) algebra.Element {
			if j >= p.settings.NumSplits() {
				return rings.NewZeroPolynomial(p.settings.BaseRing)
			}
			return p.publicParams.B.Row(j).Copy().AsVec().Dot(r).Add(sHat.Element(j))
		})
	// t[n/d] = b[n/d] * r + g
	t.SetElementAtIndex(p.settings.NumSplits(), p.publicParams.B.Row(p.settings.NumSplits()).Copy().AsVec().Dot(r).Add(g))
	// Create the masks.
	y := math.NewRandomPolynomialMatrix(p.settings.K, r.Length(), p.settings.BaseRing, p.settings.GaussianSampler)
	w := algebra.NewMatrixFromDimensions(p.settings.K, p.settings.Kappa).PopulateRows(
		func(i int) algebra.Vector {
			return p.publicParams.B0.Copy().AsMatrix().MulVec(y.Row(i))
		})
	// Create the automorphism.
	sig := math.NewAutomorphism(int64(p.settings.D), int64(p.settings.K))
	return t0, t, w, ProverState{SHat: sHat, G: g, R: r, T0: t0, Sig: sig, T: t, Y: y, W: w}
}

// CommitToRelation commits to the ternary structure of the secret and the knowledge of the secret s, s.t. As = U.
// Returns t, h, v, vp
func (p Prover) CommitToRelation(alpha algebra.Vector, gamma algebra.Matrix, state ProverState) (algebra.Vector, rings.Polynomial, rings.Polynomial, algebra.Vector, ProverState) {
	// Prover further set up.
	sum1 := CommitmentSum(p.settings.K, p.settings.NumSplits(), rings.NewPolyVec(alpha), state.Sig,
		func(i int, j int) rings.Polynomial {
			// (b[j] * y[i])^2
			tmp := p.publicParams.B.Row(j).Copy().AsVec().
				Dot(state.Y.Row(i)).
				Pow(2)
			// 3s[j] - 3
			tmp2 := state.SHat.Element(j).Copy().(rings.Polynomial).
				Scale(3).
				Sub(rings.NewOnePolynomial(p.settings.BaseRing).
					Scale(3))
			// (3s[j]-3) (b[j] * y[i])^2
			return tmp2.
				Mul(tmp).(rings.Polynomial)
		})
	sum2 := CommitmentSum(p.settings.K, p.settings.NumSplits(), rings.NewPolyVec(alpha), state.Sig,
		func(i int, j int) rings.Polynomial {
			// b[j]*y[i]
			tmp := p.publicParams.B.Row(j).Copy().AsVec().
				Dot(state.Y.Row(i))
			// 2s[j]-1
			tmp1 := state.SHat.Element(j).Copy().Scale(2).Copy().
				Sub(rings.NewOnePolynomial(p.settings.BaseRing))
			// s[j]-2
			tmp2 := state.SHat.Element(j).Copy().
				Sub(rings.NewOnePolynomial(p.settings.BaseRing).
					Scale(2))
			// (s[j]-1)s[j]
			tmp3 := state.SHat.Element(j).Copy().
				Sub(rings.NewOnePolynomial(p.settings.BaseRing)).
				Mul(state.SHat.Element(j))
			// [(2s[j]-1) (s[j]-2) + s[j](s[j]-1)](b[j]*y[i])
			return tmp.Mul(tmp1.Mul(tmp2).Add(tmp3)).(rings.Polynomial)
		})
	sum3 := CommitmentSum(p.settings.K, p.settings.NumSplits(), rings.NewPolyVec(alpha), state.Sig,
		func(i int, j int) rings.Polynomial {
			// (b[j] * y[i])^3
			return p.publicParams.B.Row(j).Copy().AsVec().
				Dot(state.Y.Row(i)).
				Pow(3).(rings.Polynomial)
		})
	state.T.SetElementAtIndex(p.settings.NumSplits()+1,
		p.publicParams.B.Row(p.settings.NumSplits()+1).Copy().AsVec().
			Dot(state.R).
			Add(p.publicParams.B.Row(p.settings.NumSplits()+2).Copy().AsVec().Dot(state.Y.Row(0))).
			Sub(sum1))
	state.T.SetElementAtIndex(p.settings.NumSplits()+2,
		p.publicParams.B.Row(p.settings.NumSplits()+2).Copy().AsVec().
			Dot(state.R).
			Add(sum2))
	v := p.publicParams.B.Row(p.settings.NumSplits() + 1).Copy().AsVec().
		Dot(state.Y.Row(0)).
		Add(sum3).(rings.Polynomial)
	At := p.publicParams.A.Copy().AsMatrix().Transpose()
	psi := algebra.NewMatrixFromDimensions(p.settings.K, p.settings.NumSplits()).PopulateRows(
		func(mu int) algebra.Vector {
			tmp := rings.NewZIntVec(At.Copy().AsMatrix().MulVec(gamma.Row(mu)))
			return SplitInvNTT(tmp, p.settings.NumSplits(), p.settings.Q, p.settings.BaseRing).AsVec()
		})
	invk := rings.NewZqInt(int64(p.settings.K), p.settings.Q).Inv().Uint64()
	gMask := LmuSum(p.settings.K, invk, state.Sig,
		func(mu int, v int) rings.Polynomial {
			// (u * gamma_mu)
			mul := p.publicParams.U.Copy().AsVec().Dot(gamma.Row(mu)).(rings.ZInt).Uint64()
			dec := rings.NewOnePolynomial(p.settings.BaseRing).Scale(mul)
			sum := algebra.NewVectorFromSize(p.settings.NumSplits()).Populate(
				func(j int) algebra.Element {
					return psi.Element(mu, j).Copy().
						Mul(state.SHat.Element(j)).(rings.Polynomial).
						Scale(uint64(p.settings.D))
				}).Sum()
			return sum.Sub(dec).(rings.Polynomial)
		})
	h := state.G.Copy().Add(gMask).(rings.Polynomial)
	vp := algebra.NewVectorFromSize(p.settings.K).Populate(
		func(i int) algebra.Element {
			// b[n/d] * y[i]
			add := p.publicParams.B.Row(p.settings.NumSplits()).Copy().AsVec().Dot(state.Y.Row(i))
			outerSum := LmuSumOuter(p.settings.K, p.settings.NumSplits(), invk, state.Sig,
				func(mu int, v int, j int) rings.Polynomial {
					index := rings.Mod(i-v, p.settings.K)
					return p.publicParams.B.Row(j).Copy().
						Mul(psi.Element(mu, j)).
						Scale(uint64(p.settings.D)).AsVec().
						Dot(state.Y.Row(index)).(rings.Polynomial)
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
func (p Prover) MaskedOpening(c rings.Polynomial, state ProverState) (algebra.Matrix, ProverState, error) {
	// Masked openings.
	z := algebra.NewMatrixFromDimensions(p.settings.K, state.R.Length()).PopulateRows(
		func(i int) algebra.Vector {
			sigc := state.Sig.Permute(int64(i), c)
			return state.Y.Row(i).Copy().Add(state.R.Copy().Mul(sigc)).AsVec()
		})

	//normInBounds := z.AllRows(
	//	func(zi math.Vector, i int) bool {
	//		infNorm := zi.NewPolyVec().InfNorm(p.settings.Q)
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
