package fastens20

import (
	"fmt"
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func extendPublicParameters(params *PublicParams, tau int) {
	extSize := tau/params.config.D + 1
	// fmt.Printf("extending public parameters with size %d\n", extSize)
	// Extend B0 horizontally by tau / d.
	b0Extension := fastmath.NewRandomPolyMatrix(params.B0.Rows(), extSize, params.config.UniformSampler, params.config.BaseRing).NTT()
	params.B0.ExtendCols(b0Extension)
	// Extend b_i vertically by tau / d.
	bExtension := fastmath.NewRandomPolyMatrix(extSize, params.B.Cols(), params.config.UniformSampler, params.config.BaseRing).NTT()
	params.B.ExtendRows(bExtension)
}

// Updates the protocol to include the approximate bound proof relation.
func updateProtocol(p *ABPProver, v *ABPVerifier, ps *ABPProverState, vs *ABPVerifierState) {
	A := p.params.A
	s := ps.S
	u := p.params.U
	// Create the new relation (A || 0, R || Id)(s, y) = (u, z)
	lrb := crypto.NewLinearRelationBuilder()
	// As = u
	A.CachedTranspose = p.params.At
	lrb.AppendEqn(crypto.NewLinearEquation(u, A.Cols()).AppendTerm(A, s))
	// Rh + y = z where h is a subvector of s

	lrb.AppendEqn(crypto.NewABPEquation(vs.ABPVerifierChal, 0, ps.ABPProverMask, ps.ABPMaskedOpening, p.params.config.BaseRing))
	newRel := lrb.Build(p.params.config.BaseRing)
	// if !newRel.IsValid() {
	// 	panic("invalid abp embedding")
	// }
	// fmt.Printf("abp equation embedded successfully: %s\n", newRel.SizesString())
	// Update the public parameters with the new relation.
	p.params.A = newRel.A
	p.params.At = newRel.At
	v.params.At = p.params.At
	v.params.A = newRel.A
	v.params.U = newRel.U
	p.params.U = newRel.U
	// Update the configuration parameters.
	p.params.config.M = newRel.A.Rows()
	v.params.config.M = newRel.A.Rows()
	// Update the N.
	p.params.config.N = newRel.A.Cols()
	v.params.config.N = newRel.A.Cols()
}

type ABPProver struct {
	Prover
	Slice fastmath.Slice
	Tau   int
}

type ABPProverState struct {
	ProverState
	ABPProverMask    *fastmath.IntVec
	ABPMaskedOpening *fastmath.IntVec
}

func NewABPProver(params PublicParams, slice fastmath.Slice, tau int) ABPProver {
	extendPublicParameters(&params, tau)
	p := NewProver(params)
	return ABPProver{
		Prover: p,
		Slice:  slice,
		Tau:    tau,
	}
}

func (p ABPProver) CommitToMessage(s *fastmath.IntVec) (*fastmath.PolyNTTVec, *fastmath.PolyNTTVec, *fastmath.PolyNTTMatrix, ABPProverState) {
	abpMask := crypto.CreateABPMask(p.Tau, p.params.config.TernarySampler, p.params.config.BaseRing)
	abpMaskPolys := abpMask.UnderlyingPolysAsPolyNTTVec()
	t0, t, w, ps := p.Prover.CommitToMessage(s)
	// Append the new commitment t_i = b_i * r + y.
	for i := 0; i < abpMaskPolys.Size(); i++ {
		t.Append(p.params.B.Row(p.params.config.NumSplits() + 3 + i).Dot(ps.R).Add(abpMaskPolys.Get(i)))
	}
	ps.T = t
	// Reperform sHat calculation.
	ps.SHat = SplitInvNTT(s.Copy().Append(abpMask), p.params).NTT()
	return t0, t, w, ABPProverState{ProverState: ps, ABPProverMask: abpMask}
}

func (p ABPProver) CreateABPMaskedOpening(abpVerifierChal *fastmath.IntMatrix, state ABPProverState) (*fastmath.IntVec, ABPProverState, error) {
	// TODO rejection sampling
	abpMaskedOpening := crypto.CreateABPMaskedOpening(abpVerifierChal, state.ABPProverMask, state.S, p.params.config.BaseRing)
	state.ABPMaskedOpening = abpMaskedOpening
	return abpMaskedOpening, state, nil
}

func (p ABPProver) CommitToRelation(alpha *fastmath.PolyNTTVec, gamma *fastmath.IntMatrix, state ABPProverState) (*fastmath.PolyNTTVec, *fastmath.PolyNTT, *fastmath.PolyNTT, *fastmath.PolyNTTVec, ABPProverState) {
	t, h, v, vp, ps := p.Prover.CommitToRelation(alpha, gamma, state.ProverState)
	state.ProverState = ps
	return t, h, v, vp, state
}

func (p ABPProver) MaskedOpening(c *fastmath.Poly, state ABPProverState) (*fastmath.PolyNTTMatrix, ABPProverState, error) {
	z, ps, err := p.Prover.MaskedOpening(c, state.ProverState)
	state.ProverState = ps
	return z, state, err
}

type ABPVerifier struct {
	Verifier
	Slice fastmath.Slice
	Tau   int
}

type ABPVerifierState struct {
	VerifierState
	ABPVerifierChal  *fastmath.IntMatrix
	ABPMaskedOpening *fastmath.IntVec
}

func NewABPVerifier(params PublicParams, slice fastmath.Slice, tau int, bound *big.Int) ABPVerifier {
	extendPublicParameters(&params, tau)
	v := NewVerifier(params)
	return ABPVerifier{
		Verifier: v,
		Slice:    slice,
		Tau:      tau,
	}
}

func (vf ABPVerifier) CreateABPChallenge() (*fastmath.IntMatrix, ABPVerifierState) {
	abpChal := crypto.CreateABPChallenge(vf.Tau, vf.params.config.N, vf.params.config.TernarySampler, vf.params.config.BaseRing)
	// Zero the irrelevant parts of R
	zeroCol := fastmath.NewIntVec(abpChal.Rows(), vf.params.config.BaseRing)
	for i := 0; i < abpChal.Cols(); i++ {
		if !vf.params.config.BoundSlice.Contains(i) {
			abpChal.SetCol(i, zeroCol)
		}
	}
	return abpChal, ABPVerifierState{ABPVerifierChal: abpChal}
}

func (vf ABPVerifier) CreateMasks(t0 *fastmath.PolyNTTVec, t *fastmath.PolyNTTVec, w *fastmath.PolyNTTMatrix, abpMaskedOpening *fastmath.IntVec, state ABPVerifierState) (*fastmath.PolyNTTVec, *fastmath.IntMatrix, ABPVerifierState) {
	alpha, gamma, vs := vf.Verifier.CreateMasks(t0, t, w)
	state.ABPMaskedOpening = abpMaskedOpening
	state.VerifierState = vs
	return alpha, gamma, state
}

func (vf ABPVerifier) CreateChallenge(t *fastmath.PolyNTTVec, h, v *fastmath.PolyNTT, vp *fastmath.PolyNTTVec, state ABPVerifierState) (*fastmath.Poly, ABPVerifierState) {
	c, vs := vf.Verifier.CreateChallenge(t, h, v, vp, state.VerifierState)
	state.VerifierState = vs
	return c, state
}

func (vf ABPVerifier) Verify(z *fastmath.PolyNTTMatrix, state ABPVerifierState) bool {
	// Infinity norm check.
	bound := vf.params.config.Bound.Uint64()
	if state.ABPMaskedOpening.Max() >= bound {
		fmt.Println("abp verifier failed infinity norm check")
		fmt.Printf("** bound = %d, norm = %d, q = %d\n", bound, state.ABPMaskedOpening.Max(), vf.params.config.Q)
		return false
	}
	res := vf.Verifier.Verify(z, state.VerifierState)
	if !res {
		fmt.Println("abp verifier original failed")
		return false
	}
	return true
}
