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
	bExtension := fastmath.NewRandomPolyMatrix(tau, params.B0.Cols(), params.config.UniformSampler, params.config.BaseRing).NTT()
	params.B.ExtendRows(bExtension)
}

// Updates the protocol to include the approximate bound proof relation.
func updateProtocol(originalRel *crypto.LinearRelation, p *ABPProver, v *ABPVerifier, ps *ABPProverState, vs *ABPVerifierState) {
	BA := originalRel.A.MulMat(vs.ABPVerifierMask)
	newRel := originalRel.Copy()
	newRel.AppendDependentOnS(BA, ps.ABPProverMask, ps.ABPMaskedOpening)
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
	ABPProverMask    *fastmath.IntVec
	ABPMaskedOpening *fastmath.IntVec
}

func NewABPProver(params PublicParams, tau int) ABPProver {
	// extendPublicParameters(&params, tau)
	p := NewProver(params)
	return ABPProver{
		Prover: p,
		Tau:    tau,
	}
}

func (p ABPProver) CommitToMessage(s *fastmath.IntVec) (*fastmath.PolyNTTVec, *fastmath.PolyNTTVec, *fastmath.PolyNTTMatrix, ABPProverState) {
	bpMask := fastmath.NewRandomTernaryIntVec(p.params.config.N, p.params.config.BaseRing)
	bpMaskPoly := fastmath.ForceNTT(bpMask.UnderlyingPolys()[0])
	t0, t, w, ps := p.Prover.CommitToMessage(s)
	for i := 1; i <= p.Tau; i++ {
		t.Append(*p.params.B.Row(p.params.config.NumSplits() + i).Dot(ps.R).Add(bpMaskPoly))
	}
	ps.SHat.Append(*bpMaskPoly)
	return t0, t, w, ABPProverState{ProverState: ps, ABPProverMask: bpMask}
}

func (p ABPProver) CreateABPMaskedOpening(bpVerifierMask *fastmath.IntMatrix, state ABPProverState) (*fastmath.IntVec, ABPProverState, error) {
	// TODO rejection sampling
	bpMaskedOpening := bpVerifierMask.MulVec(p.params.U).Add(state.ABPProverMask)
	state.ABPMaskedOpening = bpMaskedOpening
	return bpMaskedOpening, state, nil
}

func (p ABPProver) CommitToRelation(alpha *fastmath.PolyNTTVec, gamma, bpChal *fastmath.IntMatrix, state ABPProverState) (*fastmath.PolyNTTVec, *fastmath.PolyNTT, *fastmath.PolyNTT, *fastmath.PolyNTTVec, ABPProverState) {
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
	ABPVerifierMask  *fastmath.IntMatrix
	ABPMaskedOpening *fastmath.IntVec
}

func NewABPVerifier(params PublicParams, tau int) ABPVerifier {
	// extendPublicParameters(&params, tau)
	v := NewVerifier(params)
	return ABPVerifier{
		Verifier: v,
	}
}

func (vf ABPVerifier) CreateABPMask() (*fastmath.IntMatrix, ABPVerifierState) {
	bpMask := fastmath.NewRandomTernaryIntMatrix(vf.params.config.N, vf.params.config.N, vf.params.config.BaseRing)
	return bpMask, ABPVerifierState{ABPVerifierMask: bpMask}
}

func (vf ABPVerifier) CreateMasks(t0 *fastmath.PolyNTTVec, t *fastmath.PolyNTTVec, w *fastmath.PolyNTTMatrix, bpMaskedOpening *fastmath.IntVec, state ABPVerifierState) (*fastmath.PolyNTTVec, *fastmath.IntMatrix, ABPVerifierState) {
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
	if state.ABPMaskedOpening.InfNorm() >= bound {
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
	bpVerifierMask, vs := verifier.CreateABPMask()
	bpMaskedOpening, ps, _ := prover.CreateABPMaskedOpening(bpVerifierMask, ps)
	// Update the relation to embed the approximate bound proof.
	updateProtocol(&linRel, &prover, &verifier, &ps, &vs)
	// Resume the normal execution.
	alpha, gamma, vs := verifier.CreateMasks(t0, t, w, bpMaskedOpening, vs)
	// Continue the execution.
	t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, bpVerifierMask, ps)
	c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
	z, ps, _ := prover.MaskedOpening(c, ps)
	return verifier.Verify(z, vs)
}
