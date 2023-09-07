package fastens20

import (
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
)

type ProverState struct {
	G    *fastmath.PolyNTT
	S    *fastmath.IntVec
	SHat *fastmath.PolyNTTVec
	R    *fastmath.PolyNTTVec
	T0   *fastmath.PolyNTTVec
	T    *fastmath.PolyNTTVec
	Y    *fastmath.PolyNTTMatrix
	W    *fastmath.PolyNTTMatrix
	V    *fastmath.PolyNTT
	Psi  *fastmath.PolyNTTMatrix
	H    *fastmath.PolyNTT
	Vp   *fastmath.PolyNTTVec
	Z    *fastmath.PolyNTTMatrix
}

type Prover struct {
	params PublicParams
}

func NewProver(params PublicParams) Prover {
	return Prover{params}
}

// CommitToMessage commits to the given secret message s.
// Returns t0, t, w
func (p Prover) CommitToMessage(s *fastmath.IntVec) (*fastmath.PolyNTTVec, *fastmath.PolyNTTVec, *fastmath.PolyNTTMatrix, ProverState) {
	e := logging.LogExecStart("Prover.CommitToMessage", "working")
	defer e.LogExecEnd()
	// Rebase the message into polynomial space.
	sHat := SplitInvNTT(s, p.params).NTT()
	// Sample a polynomial g s.t. g_0=...=g_{k-1}=0
	gPoly := fastmath.NewRandomPoly(p.params.config.UniformSampler, p.params.config.BaseRing)
	for i := 0; i < p.params.config.K; i++ {
		gPoly.SetForce(i, 0)
	}
	g := gPoly.NTT()
	// Sample the randomness.
	r := fastmath.NewRandomPolyVec(p.params.B0.Cols(), p.params.config.TernarySampler, p.params.config.BaseRing).NTT()
	// Compute the commitments.
	t0 := p.params.B0.MulVec(r)
	t := fastmath.NewPolyVec(p.params.config.NumSplits(), p.params.config.BaseRing).NTT()
	t.Populate(func(j int) *fastmath.PolyNTT {
		return p.params.B.Row(j).Dot(r).Add(sHat.Get(j))
	})
	// t[n/d + 1] = b[n/d + 1] * r + g
	t.Append(p.params.B.Row(p.params.config.NumSplits()).Dot(r).Add(g))
	// Create the masks.
	y := fastmath.NewRandomPolyMatrix(p.params.config.K, r.Size(), p.params.config.GaussianSampler, p.params.config.BaseRing).NTT()
	w := fastmath.NewPolyMatrix(p.params.config.K, p.params.config.Kappa, p.params.config.BaseRing).NTT()
	w.PopulateRows(func(i int) *fastmath.PolyNTTVec {
		return p.params.B0.MulVec(y.Row(i))
	})
	return t0.Copy(), t.Copy(), w.Copy(), ProverState{S: s.Copy(), SHat: sHat, G: g, R: r, T0: t0, T: t, Y: y, W: w}
}

// CommitToRelation commits to the ternary structure of the secret and the knowledge of the secret s, s.t. As = U.
// Returns t, h, v, vp
func (p Prover) CommitToRelation(alpha *fastmath.PolyNTTVec, gamma *fastmath.IntMatrix, state ProverState) (*fastmath.PolyNTTVec, *fastmath.PolyNTT, *fastmath.PolyNTT, *fastmath.PolyNTTVec, ProverState) {
	e0 := logging.LogExecStart("Prover.CommitToRelation", "working")
	defer e0.LogExecEnd()
	// Get the ternary part of the s.
	e := logging.LogExecStart("Prover.CommitToRelation", "commitment sums")
	sum1 := CommitmentSum(p.params.config.K, p.params.config.NumTernarySplits(), alpha,
		func(i int, j int) *fastmath.PolyNTT {
			jp := p.params.config.TernarySlice.Start/p.params.config.D + j
			// (b[j]*y[i])^2
			tmp := p.params.B.Row(jp).Dot(state.Y.Row(i)).
				Pow(2)
			// 3s[j]-3
			neg1 := fastmath.NewOnePoly(1, p.params.config.BaseRing).NTT().Neg()
			tmp2 := neg1.Add(state.SHat.Get(jp)).Scale(3)
			// (3s[j]-3) * (b[j]*y[i])^2
			return tmp2.Mul(tmp)
		}, p.params)
	sum2 := CommitmentSum(p.params.config.K, p.params.config.NumTernarySplits(), alpha,
		func(i int, j int) *fastmath.PolyNTT {
			jp := p.params.config.TernarySlice.Start/p.params.config.D + j
			// b[j]*y[i]
			tmp := p.params.B.Row(jp).Dot(state.Y.Row(i))
			// 2s[j]-1
			neg1 := fastmath.NewOnePoly(1, p.params.config.BaseRing).NTT().Neg()
			tmp1 := state.SHat.Get(jp).Copy().Scale(2).Add(neg1)
			// s[j]-2
			tmp2 := fastmath.NewOnePoly(2, p.params.config.BaseRing).NTT().Neg().Add(state.SHat.Get(jp))
			// (s[j]-1)*s[j]
			tmp3 := neg1.Add(state.SHat.Get(jp)).Mul(state.SHat.Get(jp))
			// [(2s[j]-1) * (s[j]-2) + s[j](s[j]-1)] * (b[j]*y[i])
			return tmp.Mul(tmp1.Mul(tmp2).Add(tmp3))
		}, p.params)
	sum3 := CommitmentSum(p.params.config.K, p.params.config.NumTernarySplits(), alpha,
		func(i int, j int) *fastmath.PolyNTT {
			jp := p.params.config.TernarySlice.Start/p.params.config.D + j
			// (b[j] * y[i])^3
			by := p.params.B.Row(jp).Dot(state.Y.Row(i))
			return by.Pow(3)
		}, p.params)
	// t[n/d+2]
	state.T.Append(
		p.params.B.Row(p.params.config.NumSplits() + 1).
			Dot(state.R).
			Add(p.params.B.Row(p.params.config.NumSplits() + 2).
				Dot(state.Y.Row(0))).
			Add(sum1.Neg()))
	// t[n/d+3]
	state.T.Append(
		p.params.B.Row(p.params.config.NumSplits() + 2).
			Dot(state.R).
			Add(sum2))
	e.LogExecEnd()

	e = logging.LogExecStart("Prover.CommitToRelation", "v calculation")
	v := p.params.B.Row(p.params.config.NumSplits() + 1).Dot(state.Y.Row(0)).Add(sum3)
	e.LogExecEnd()

	e = logging.LogExecStart("Prover.CommitToRelation", "psi calculation")
	psi := fastmath.NewPolyMatrix(p.params.config.K, p.params.config.NumSplits(), p.params.config.BaseRing).NTT()
	psi.PopulateRows(func(mu int) *fastmath.PolyNTTVec {
		tmp := p.params.A.MulVecTranspose(gamma.RowView(mu))
		return SplitInvNTT(tmp, p.params).NTT()
	})
	e.LogExecEnd()

	e = logging.LogExecStart("Prover.CommitToRelation", "h calculation")
	gMask := LmuSum(p.params.config.K, p.params.config.Cache.InvK,
		func(mu int, v int) *fastmath.PolyNTT {
			presum := fastmath.NewPolyVec(p.params.config.NumSplits(), p.params.config.BaseRing).NTT()
			presum.Populate(func(j int) *fastmath.PolyNTT {
				return psi.Get(mu, j).Copy().
					Mul(state.SHat.Get(j)).
					Scale(uint64(p.params.config.D))
			})
			// (u * gamma_mu)
			dec := fastmath.NewOnePolyLevels(p.params.U.Dot(gamma.RowView(mu)), p.params.config.BaseRing).NTT()
			return presum.Sum().Add(dec.Neg())
		}, p.params)
	h := gMask.Add(state.G)
	e.LogExecEnd()

	e = logging.LogExecStart("Prover.CommitToRelation", "vp calculation")
	vp := fastmath.NewPolyVec(p.params.config.K, p.params.config.BaseRing).NTT()
	vp.Populate(func(i int) *fastmath.PolyNTT {
		// b[n/d + 1] * y[i]
		add := p.params.B.Row(p.params.config.NumSplits()).Dot(state.Y.Row(i))
		outerSum := LmuSumOuter(p.params.config.K, p.params.config.NumSplits(), p.params.config.Cache.InvK,
			func(mu int, v int, j int) *fastmath.PolyNTT {
				index := i - v
				if index < 0 {
					index += p.params.config.K
				}
				index = index % p.params.config.K
				dotResult := p.params.B.Row(j).Dot(state.Y.Row(index))
				return dotResult.Scale(uint64(p.params.config.D)).Mul(psi.Get(mu, j))
			}, p.params)
		return outerSum.Add(add)
	})
	e.LogExecEnd()
	// Update the state.
	state.V = v
	state.Psi = psi
	state.H = h
	state.Vp = vp
	return state.T.Copy(), state.H.Copy(), v.Copy(), vp.Copy(), state
}

// MaskedOpening returns the masked openings to the commitments.
// Returns z
func (p Prover) MaskedOpening(c *fastmath.Poly, state ProverState) (*fastmath.PolyNTTMatrix, ProverState, error) {
	e := logging.LogExecStart("Prover.MaskedOpening", "working")
	defer e.LogExecEnd()
	// Masked openings.
	z := fastmath.NewPolyMatrix(p.params.config.K, state.R.Size(), p.params.config.BaseRing).NTT()
	z.PopulateRows(func(i int) *fastmath.PolyNTTVec {
		sigc := c.Permute(int64(i), p.params.Sig).NTT()
		tmp := state.R.Copy().MulAll(sigc)
		return tmp.Add(state.Y.Row(i))
	})
	state.Z = z
	return z.Copy(), state, nil
}
