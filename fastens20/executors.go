package fastens20

import (
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
)

func executeWithoutBoundProof(s *fastmath.IntVec, params PublicParams) bool {
	prover := NewProver(params)
	verifier := NewVerifier(params)
	// Commit to the message.
	e := logging.LogExecStart("Execute", "Prover.CommitToMessage")
	t0, t, w, ps := prover.CommitToMessage(s)
	e.LogExecEnd()

	e = logging.LogExecStart("Execute", "Verifier.CreateMasks")
	alpha, gamma, vs := verifier.CreateMasks(t0, t, w)
	e.LogExecEnd()

	e = logging.LogExecStart("Execute", "Prover.CommitToRelation")
	t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
	e.LogExecEnd()

	e = logging.LogExecStart("Execute", "Verifier.CreateChallenge")
	c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
	e.LogExecEnd()

	// Recreate the masked opening until it satisfies the shortness condition.
	e = logging.LogExecStart("Execute", "Prover.MaskedOpening")
	z, ps, err := prover.MaskedOpening(c, ps)
	for err != nil {
		z, ps, err = prover.MaskedOpening(c, ps)
	}
	e.LogExecEnd()

	e = logging.LogExecStart("Execute", "Verifier.Verify")
	res := verifier.Verify(z, vs)
	e.LogExecEnd()

	return res
}

func executeWithBoundProof(s *fastmath.IntVec, params PublicParams) bool {
	prover := NewABPProver(params, params.config.BoundSlice, params.config.Tau)
	verifier := NewABPVerifier(params, params.config.BoundSlice, params.config.Tau, params.config.Bound)
	// fmt.Println("abp exchange initiated")
	// Commit to the message.
	e := logging.LogExecStart("Execute", "Prover.CommitToMessage")
	t0, t, w, ps := prover.CommitToMessage(s)
	e.LogExecEnd()
	// ABP exchange.
	e = logging.LogExecStart("Execute", "Verifier.CreateABPChallenge")
	abpVerifierChal, vs := verifier.CreateABPChallenge()
	e.LogExecEnd()

	e = logging.LogExecStart("Execute", "Prover.CreateABPMaskedOpening")
	abpMaskedOpening, ps, _ := prover.CreateABPMaskedOpening(abpVerifierChal, ps)
	e.LogExecEnd()
	// fmt.Println("updating the protocol")
	// Update the relation to embed the approximate bound proof.
	logging.LogExecution("Execute", "updating", func() interface{} { updateProtocol(&prover, &verifier, &ps, &vs); return nil })
	// fmt.Println("continuing the executing of the protocol")
	// Resume the normal execution.
	e = logging.LogExecStart("Execute", "Verifier.CreateMasks")
	alpha, gamma, vs := verifier.CreateMasks(t0, t, w, abpMaskedOpening, vs)
	e.LogExecEnd()
	// Continue the execution.
	e = logging.LogExecStart("Execute", "Prover.CommitToRelation")
	t, h, v, vp, ps := prover.CommitToRelation(alpha, gamma, ps)
	e.LogExecEnd()

	e = logging.LogExecStart("Execute", "Verifier.CreateChallenge")
	c, vs := verifier.CreateChallenge(t, h, v, vp, vs)
	e.LogExecEnd()

	e = logging.LogExecStart("Execute", "Prover.MaskedOpening")
	z, ps, _ := prover.MaskedOpening(c, ps)
	e.LogExecEnd()

	e = logging.LogExecStart("Execute", "Verifier.Verify")
	res := verifier.Verify(z, vs)
	e.LogExecEnd()

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