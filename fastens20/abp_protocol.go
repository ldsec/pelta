package fastens20

import (
	"fmt"
	"math/big"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func extendPublicParameters(params *PublicParams, tau int) {
	b0Extension := fastmath.NewRandomPolyMatrix(params.config.Kappa, tau, params.config.UniformSampler, params.config.BaseRing).NTT()
	params.B0.ExtendCols(b0Extension)
	bExtension := fastmath.NewRandomPolyMatrix(tau, params.B.Cols(), params.config.UniformSampler, params.config.BaseRing).NTT()
	params.B.ExtendRows(bExtension)
}

// Updates the protocol to include the approximate bound proof relation.
func updateProtocol(originalRel *crypto.LinearRelation, p *ABPProver, v *ABPVerifier, ps *ABPProverState, vs *ABPVerifierState) {
	BA := vs.ABPVerifierChal.ToIntMatrix().MulMat(originalRel.A)
	y := ps.ABPProverMask.ToIntVec()
	z := ps.ABPMaskedOpening.ToIntVec()
	newRel := originalRel.Copy()
	newRel.AppendDependentOnS(BA, y, z)
	p.params.A = newRel.A
	v.params.A = newRel.A
	v.params.U = newRel.U
	p.params.U = newRel.U
	p.params.config.M = newRel.A.Rows()
	v.params.config.M = newRel.A.Rows()
	p.params.config.N = newRel.A.Cols()
	v.params.config.N = newRel.A.Cols()
}

type ABPProver struct {
	Prover
	Tau int
}

type ABPProverState struct {
	ProverState
	ABPProverMask    *fastmath.PolyNTTVec
	ABPMaskedOpening *fastmath.PolyNTTVec
}

func NewABPProver(params PublicParams, tau int) ABPProver {
	extendPublicParameters(&params, tau)
	p := NewProver(params)
	return ABPProver{
		Prover: p,
		Tau:    tau,
	}
}

func (p ABPProver) CommitToMessage(s *fastmath.IntVec) (*fastmath.PolyNTTVec, *fastmath.PolyNTTVec, *fastmath.PolyNTTMatrix, ABPProverState) {
	abpMask := fastmath.NewRandomTernaryPolyVec(p.Tau, p.params.config.BaseRing)
	abpMaskNTT := fastmath.NewPolyVec(p.Tau, p.params.config.BaseRing).NTT()
	for i := 0; i < p.Tau; i++ {
		abpMaskNTT.Set(i, fastmath.ForceNTT(abpMask.Get(i)))
	}
	t0, t, w, ps := p.Prover.CommitToMessage(s)
	for i := 0; i < p.Tau; i++ {
		t.Append(p.params.B.Row(p.params.config.NumSplits() + 1 + i).Dot(ps.R).Add(abpMaskNTT.Get(i)))
		ps.SHat.Append(abpMaskNTT.Get(i))
	}
	return t0, t, w, ABPProverState{ProverState: ps, ABPProverMask: abpMaskNTT}
}

func (p ABPProver) CreateABPMaskedOpening(bpVerifierChal *fastmath.PolyNTTMatrix, state ABPProverState) (*fastmath.PolyNTTVec, ABPProverState, error) {
	// TODO rejection sampling
	abpMaskedOpening := bpVerifierChal.ToIntMatrix().MulVec(p.params.U).Add(state.ABPProverMask.ToIntVec()).UnderlyingPolysAsPolyNTTVec()
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
	Tau int
}

type ABPVerifierState struct {
	VerifierState
	ABPVerifierChal  *fastmath.PolyNTTMatrix
	ABPMaskedOpening *fastmath.PolyNTTVec
}

func NewABPVerifier(params PublicParams, tau int) ABPVerifier {
	extendPublicParameters(&params, tau)
	v := NewVerifier(params)
	return ABPVerifier{
		Verifier: v,
		Tau:      tau,
	}
}

func (vf ABPVerifier) CreateABPChallenge() (*fastmath.PolyNTTMatrix, ABPVerifierState) {
	l := vf.params.U.Size() / vf.params.config.D
	abpChal := fastmath.NewRandomPolyMatrix(vf.Tau*vf.params.config.D, l, vf.params.config.TernarySampler, vf.params.config.BaseRing).NTT()
	return abpChal, ABPVerifierState{ABPVerifierChal: abpChal}
}

func (vf ABPVerifier) CreateMasks(t0 *fastmath.PolyNTTVec, t *fastmath.PolyNTTVec, w *fastmath.PolyNTTMatrix, bpMaskedOpening *fastmath.PolyNTTVec, state ABPVerifierState) (*fastmath.PolyNTTVec, *fastmath.IntMatrix, ABPVerifierState) {
	alpha, gamma, vs := vf.Verifier.CreateMasks(t0, t, w)
	state.ABPMaskedOpening = bpMaskedOpening
	state.VerifierState = vs
	return alpha, gamma, state
}

func (vf ABPVerifier) CreateChallenge(t *fastmath.PolyNTTVec, h, v *fastmath.PolyNTT, vp *fastmath.PolyNTTVec, state ABPVerifierState) (*fastmath.Poly, ABPVerifierState) {
	c, vs := vf.Verifier.CreateChallenge(t, h, v, vp, state.VerifierState)
	state.VerifierState = vs
	return c, state
}

func (vf ABPVerifier) Verify(z *fastmath.PolyNTTMatrix, state ABPVerifierState) bool {
	bound := big.NewInt(0).Sqrt(
		big.NewInt(0).Div(
			vf.params.config.Q,
			big.NewInt(int64(2*vf.params.U.Size())))).Uint64() / 2
	if state.ABPMaskedOpening.Copy().InvNTT().InfNorm() >= bound {
		fmt.Println("verifier failed ABP check")
		return false
	}
	return vf.Verifier.Verify(z, state.VerifierState)
}

func ExecuteWithBoundProof(s *fastmath.IntVec, tau int, params PublicParams) bool {
	linRel := crypto.NewLinearRelation(params.A, s)
	prover := NewABPProver(params, tau)
	verifier := NewABPVerifier(params, tau)
	// Commit to the message.
	t0, t, w, ps := prover.CommitToMessage(s)
	// ABP exchange.
	bpVerifierChal, vs := verifier.CreateABPChallenge()
	bpMaskedOpening, ps, _ := prover.CreateABPMaskedOpening(bpVerifierChal, ps)
	// Update the relation to embed the approximate bound proof.
	updateProtocol(&linRel, &prover, &verifier, &ps, &vs)
	// Resume the normal execution.
	alpha, gamma, vs := verifier.CreateMasks(t0, t, w, bpMaskedOpening, vs)
	// Continue the execution.
	t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
	c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
	z, ps, _ := prover.MaskedOpening(c, ps)
	return verifier.Verify(z, vs)
}
