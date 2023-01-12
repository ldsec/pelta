package fastens20

import (
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
)

func executeWithoutBoundProof(s *fastmath.IntVec, params PublicParams) bool {
	prover := NewProver(params)
	verifier := NewVerifier(params)
	// Commit to the message.
	t0, t, w, ps := prover.CommitToMessage(s)
	alpha, gamma, vs := verifier.CreateMasks(t0, t, w)
	t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
	c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
	// Recreate the masked opening until it satisfies the shortness condition.
	z, ps, err := prover.MaskedOpening(c, ps)
	for err != nil {
		z, ps, err = prover.MaskedOpening(c, ps)
	}
	res := verifier.Verify(z, vs)
	return res
}

func executeWithBoundProof(s *fastmath.IntVec, params PublicParams) bool {
	prover := NewABPProver(params, params.config.BoundSlice, params.config.Tau)
	verifier := NewABPVerifier(params, params.config.BoundSlice, params.config.Tau, params.config.Bound)
	// fmt.Println("abp exchange initiated")
	// Commit to the message.
	t0, t, w, ps := prover.CommitToMessage(s)
	// ABP exchange.
	abpVerifierChal, vs := verifier.CreateABPChallenge()
	abpMaskedOpening, ps, _ := prover.CreateABPMaskedOpening(abpVerifierChal, ps)
	// fmt.Println("updating the protocol")
	// Update the relation to embed the approximate bound proof.
	updateProtocol(&prover, &verifier, &ps, &vs)
	// fmt.Println("continuing the executing of the protocol")
	// Resume the normal execution.
	alpha, gamma, vs := verifier.CreateMasks(t0, t, w, abpMaskedOpening, vs)
	// Continue the execution.
	t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
	c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
	z, ps, _ := prover.MaskedOpening(c, ps)
	res := verifier.Verify(z, vs)
	return res
}

// Execute runs the augmented ENS20 protocol.
func Execute(s *fastmath.IntVec, params PublicParams) bool {
	// Either execute with or without bound proof.
	var res bool
	e := logging.LogExecStart("Execute", "Protocol execution")
	if params.config.ABPEnabled {
		res = executeWithBoundProof(s, params)
	} else {
		res = executeWithoutBoundProof(s, params)
	}
	e.LogExecEnd()
	return res
}