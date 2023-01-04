package fastens20

import (
	"fmt"
	"time"

	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

func ExecuteWithoutBoundProof(s *fastmath.IntVec, params PublicParams) bool {
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

func ExecuteWithBoundProof(s *fastmath.IntVec, slice fastmath.Slice, params PublicParams) bool {
	prover := NewABPProver(params, slice, params.config.Tau)
	verifier := NewABPVerifier(params, slice, params.config.Tau)
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
	return verifier.Verify(z, vs)
}

// Execute runs the augmented ENS20 protocol. `ter` and `abp` proof definitions are optional.
func Execute(s *fastmath.IntVec, val *crypto.ValidityProofDef, ter *crypto.TernaryProofDef, abp *crypto.ApproxBoundProofDef, ringParams fastmath.RingParams) bool {
	// Generate the public parameters
	config := DefaultProtocolConfig(ringParams, val.Rel).WithTernarySlice(fastmath.NewSlice(0, 0))
	if ter != nil {
		config = config.WithTernarySlice(ter.Target)
	}
	params := GeneratePublicParameters(config, val.Rel)
	// Either execute with or without bound proof.
	var res bool
	t0 := time.Now()
	if abp != nil {
		res = ExecuteWithBoundProof(s, abp.Target, params)
	} else {
		res = ExecuteWithoutBoundProof(s, params)
	}
	dt := time.Now().Sub(t0)
	fmt.Printf("protocol execution took %dms\n", dt.Milliseconds())
	return res
}