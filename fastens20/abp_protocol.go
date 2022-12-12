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
func updateProtocol(p *ABPProver, v *ABPVerifier, ps *ABPProverState, vs *ABPVerifierState) {
	A := p.params.A
	s := ps.S
	u := p.params.U
	lrb := crypto.NewLinearRelationBuilder()
	lrb.AppendEqn(crypto.NewLinearEquation(u, A.Cols()).AppendTerm(A, s))
	lrb.AppendEqn(crypto.NewABPEquation(vs.ABPVerifierChal, A, 0, ps.ABPProverMask, ps.ABPMaskedOpening, p.params.config.BaseRing))
	fmt.Println("building")
	newRel := lrb.Build(p.params.config.BaseRing)
	if !newRel.Verify() {
		fmt.Println("invalid ABP embedding")
	}
	p.params.A = newRel.A
	v.params.A = newRel.A
	v.params.U = newRel.U
	p.params.U = newRel.U
	p.params.config.M = newRel.A.Rows()
	v.params.config.M = newRel.A.Rows()
	p.params.config.TernaryLength = p.params.config.N
	v.params.config.TernaryLength = v.params.config.N
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
	abpMask := crypto.CreateABPMask(p.Tau, p.params.config.TernarySampler, p.params.config.BaseRing)
	t0, t, w, ps := p.Prover.CommitToMessage(s)
	for i := 0; i < p.Tau; i++ {
		t.Append(p.params.B.Row(p.params.config.NumSplits() + 1 + i).Dot(ps.R).Add(abpMask.Get(i)))
		ps.SHat.Append(abpMask.Get(i))
	}
	ps.T = t
	return t0, t, w, ABPProverState{ProverState: ps, ABPProverMask: abpMask}
}

func (p ABPProver) CreateABPMaskedOpening(bpVerifierChal *fastmath.PolyNTTMatrix, state ABPProverState) (*fastmath.PolyNTTVec, ABPProverState, error) {
	// TODO rejection sampling
	abpMaskedOpening := crypto.CreateABPMaskedOpening(bpVerifierChal, state.ABPProverMask, p.params.U, p.params.config.BaseRing)
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
	abpChal := crypto.CreateABPChallenge(vf.Tau, vf.params.config.M, vf.params.config.TernarySampler, vf.params.config.BaseRing)
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
	res := vf.Verifier.Verify(z, state.VerifierState)
	if !res {
		fmt.Println("verifier original failed")
		return false
	}
	bound := big.NewInt(0).Sqrt(
		big.NewInt(0).Div(
			vf.params.config.Q,
			big.NewInt(int64(2*vf.params.U.Size())))).Uint64() / 2
	if state.ABPMaskedOpening.Copy().InvNTT().InfNorm() >= bound {
		fmt.Println("verifier failed ABP check")
		return false
	}
	return true
}

func ExecuteWithBoundProof(s *fastmath.IntVec, tau int, params PublicParams) bool {
	prover := NewABPProver(params, tau)
	verifier := NewABPVerifier(params, tau)
	// Commit to the message.
	t0, t, w, ps := prover.CommitToMessage(s)
	// ABP exchange.
	bpVerifierChal, vs := verifier.CreateABPChallenge()
	bpMaskedOpening, ps, _ := prover.CreateABPMaskedOpening(bpVerifierChal, ps)
	fmt.Println("executing the protocol (1)")
	// Update the relation to embed the approximate bound proof.
	updateProtocol(&prover, &verifier, &ps, &vs)
	fmt.Println("executing the protocol (2)")
	// Resume the normal execution.
	alpha, gamma, vs := verifier.CreateMasks(t0, t, w, bpMaskedOpening, vs)
	// Continue the execution.
	t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
	c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
	z, ps, _ := prover.MaskedOpening(c, ps)
	return verifier.Verify(z, vs)
}
