package fastens20

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

type ProverState struct {
	G    fastmath.Poly
	SHat fastmath.PolyVec
	Sig  fastmath.Automorphism
	R    fastmath.PolyVec
	T0   fastmath.PolyVec
	T    fastmath.PolyVec
	Y    fastmath.PolyMatrix
	W    fastmath.PolyMatrix
	V    fastmath.Poly
	Psi  fastmath.PolyMatrix
	H    fastmath.Poly
	Vp   fastmath.PolyVec
	Z    fastmath.PolyMatrix
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
func (p Prover) CommitToMessage(s fastmath.IntVec) (fastmath.PolyVec, fastmath.PolyVec, fastmath.PolyMatrix, ProverState) {
	// Split the message into polynomial space.
	sHat := SplitInvNTT(s, p.settings.BaseRing)
	// Sample a polynomial g s.t. g_0=...=g_{K-1}=0
	g := fastmath.NewRandomPoly(p.settings.UniformSampler, p.settings.BaseRing)
	for i := 0; i < p.settings.K; i++ {
		g.Set(i, 0)
	}
	// Sample the randomness.
	r := fastmath.NewRandomPolyVec(p.publicParams.B0.Cols(), p.settings.TernarySampler, p.settings.BaseRing)
	// Compute the commitments.
	br := fastmath.NewPolyVec(p.publicParams.B.Rows(), p.settings.BaseRing)
	br.Populate(func(i int) fastmath.Poly {
		return *p.publicParams.B.Row(i).Dot(&r)
	})
	t0 := p.publicParams.B0.MulVec(&r)
	t := fastmath.NewPolyVec(p.settings.NumSplits(), p.settings.BaseRing)
	t.Populate(func(j int) fastmath.Poly {
		return *br.Get(j).Add(sHat.Get(j))
	})
	// t[n/d] = b[n/d] * r + g
	t.Append(*br.Get(p.settings.NumSplits()).Add(&g))
	// Create the masks.
	y := fastmath.NewRandomPolyMatrix(p.settings.K, r.Size(), p.settings.GaussianSampler, p.settings.BaseRing)
	w := fastmath.NewPolyMatrix(p.settings.K, p.settings.Kappa, p.settings.BaseRing)
	w.PopulateRows(func(i int) fastmath.PolyVec {
		return p.publicParams.B0.MulVec(y.Row(i))
	})
	// Create the automorphism.
	sig := fastmath.NewAutomorphism(uint64(p.settings.D), uint64(p.settings.K))
	return t0, t, w, ProverState{SHat: sHat, G: g, R: r, T0: t0, Sig: sig, T: t, Y: y, W: w}
}

// CommitToRelation commits to the ternary structure of the secret and the knowledge of the secret s, s.t. As = U.
// Returns t, h, v, vp
func (p Prover) CommitToRelation(alpha fastmath.PolyVec, gamma fastmath.IntMatrix, state ProverState) (fastmath.PolyVec, fastmath.Poly, fastmath.Poly, fastmath.PolyVec, ProverState) {
	// Prover further set up.
	sum1 := CommitmentSum(p.settings.K, p.settings.NumSplits(), alpha, state.Sig,
		func(i int, j int) fastmath.Poly {
			// (b[j] * y[i])^2
			tmp := p.publicParams.B.Row(j).Dot(state.Y.Row(i))
			tmp.PowModCoeffs(2, p.settings.Q.Uint64())
			// 3s[j] - 3
			tmp2 := state.SHat.Get(j).Copy()
			three := fastmath.NewOnePoly(3, p.settings.BaseRing)
			tmp2.Scale(3).Add(three.Neg())
			// (3s[j]-3) (b[j] * y[i])^2
			return *tmp2.MulCoeffs(tmp)
		}, p.settings.BaseRing)
	sum2 := CommitmentSum(p.settings.K, p.settings.NumSplits(), alpha, state.Sig,
		func(i int, j int) fastmath.Poly {
			// b[j]*y[i]
			tmp := p.publicParams.B.Row(j).Dot(state.Y.Row(i))
			// 2s[j]-1
			tmp1 := state.SHat.Get(j).Copy()
			one := fastmath.NewOnePoly(1, p.settings.BaseRing)
			tmp1.Scale(2).Add(one.Neg())
			// s[j]-2
			tmp2 := state.SHat.Get(j).Copy()
			two := fastmath.NewOnePoly(2, p.settings.BaseRing)
			tmp2.Add(two.Neg())
			// (s[j]-1)s[j]
			tmp3 := state.SHat.Get(j).Copy().Add(&one).MulCoeffs(state.SHat.Get(j))
			// [(2s[j]-1) (s[j]-2) + s[j](s[j]-1)](b[j]*y[i])
			return *tmp.MulCoeffs(tmp1.MulCoeffs(tmp2).Add(tmp3))
		}, p.settings.BaseRing)
	sum3 := CommitmentSum(p.settings.K, p.settings.NumSplits(), alpha, state.Sig,
		func(i int, j int) fastmath.Poly {
			// (b[j] * y[i])^3
			by := p.publicParams.B.Row(j).Dot(state.Y.Row(i))
			return *by.PowModCoeffs(3, p.settings.Q.Uint64())
		}, p.settings.BaseRing)
	state.T.Append(*p.publicParams.B.Row(p.settings.NumSplits() + 1).
		Dot(&state.R).
		Add(p.publicParams.B.Row(p.settings.NumSplits() + 2).
			Dot(state.Y.Row(0))).
		Add(sum1.Neg()))
	state.T.Append(*p.publicParams.B.Row(p.settings.NumSplits() + 2).
		Dot(&state.R).
		Add(&sum2))
	v := *p.publicParams.B.Row(p.settings.NumSplits() + 1).Dot(state.Y.Row(0)).Add(&sum3)
	At := p.publicParams.A.Transposed()
	psi := fastmath.NewPolyMatrix(p.settings.K, p.settings.NumSplits(), p.settings.BaseRing)
	psi.PopulateRows(func(mu int) fastmath.PolyVec {
		tmp := At.MulVec(gamma.RowView(mu))
		return SplitInvNTT(tmp, p.settings.BaseRing)
	})
	invk := big.NewInt(0).ModInverse(big.NewInt(int64(p.settings.K)), p.settings.Q).Uint64()
	gMask := LmuSum(p.settings.K, invk, state.Sig,
		func(mu int, v int) fastmath.Poly {
			// (u * gamma_mu)
			mul := p.publicParams.U.Dot(gamma.RowView(mu))
			dec := fastmath.NewOnePoly(mul, p.settings.BaseRing)
			presum := fastmath.NewPolyVec(p.settings.NumSplits(), p.settings.BaseRing)
			presum.Populate(func(j int) fastmath.Poly {
				return *psi.Get(mu, j).
					MulCoeffs(state.SHat.Get(j)).
					Scale(uint64(p.settings.D))
			})
			return *presum.Sum().Add(dec.Neg())
		}, p.settings.BaseRing)
	h := state.G.Copy().Add(&gMask)
	vp := fastmath.NewPolyVec(p.settings.K, p.settings.BaseRing)
	vp.Populate(func(i int) fastmath.Poly {
		// b[n/d] * y[i]
		add := p.publicParams.B.Row(p.settings.NumSplits()).Dot(state.Y.Row(i))
		outerSum := LmuSumOuter(p.settings.K, p.settings.NumSplits(), invk, state.Sig,
			func(mu int, v int, j int) fastmath.Poly {
				index := big.NewInt(0).
					Mod(big.NewInt(int64(i-v)),
						big.NewInt(int64(p.settings.K))).
					Uint64()
				return *p.publicParams.B.Row(j).
					MulAll(psi.Get(mu, j)).
					ScaleAll(uint64(p.settings.D)).
					Dot(state.Y.Row(int(index)))
			}, p.settings.BaseRing)
		return *outerSum.Add(add)
	})
	// Update the state.
	state.V = v
	state.Psi = psi
	state.H = *h
	state.Vp = vp
	return state.T, state.H, v, vp, state
}

// MaskedOpening returns the masked openings to the commitments.
// Returns z
func (p Prover) MaskedOpening(c fastmath.Poly, state ProverState) (fastmath.PolyMatrix, ProverState, error) {
	// Masked openings.
	z := fastmath.NewPolyMatrix(p.settings.K, state.R.Size(), p.settings.BaseRing)
	z.PopulateRows(func(i int) fastmath.PolyVec {
		sigc := c.Permute(int64(i), state.Sig)
		tmp := state.R.Copy().MulAll(&sigc)
		return *state.Y.Row(i).Copy().Add(tmp)
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
