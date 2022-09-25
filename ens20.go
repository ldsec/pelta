package main

import "github.com/ldsec/lattigo/v2/ring"

// Prover represents the prover party in the protocol.
type Prover struct{}

// Verifier represents the verifier party in the protocol.
type Verifier struct{}

// Protocol represents an instance of the protocol described in ENS20.
// ...
type Protocol struct {
	Prover   Prover
	Verifier Verifier
}

type PublicParams struct {
	N         uint64 // dimension of Rq
	n         uint64 // dimension for MSIS
	m1        uint64 // dimension of the input in Rq
	m2        uint64 // dimension of the input in Rq
	lbd       uint64 // dimension for MLWE
	sizeB     uint64 // dimension of bi
	size_t    uint64 // dimension of the commitment
	size_bVec uint64 // dimension of bVec vector of {bi}
	k         uint64 // repeat rate in ENS20
	d1        uint64 // delta1 bound on the mask // Change for appropriate value
}

// Returns the commitments t, w_i.
func (*Prover) initiate() ([]*ring.Poly, [][]*ring.Poly, error) {
	return nil, nil, nil
}

// Returns the gamma values.
func (*Verifier) createGamma() ([][]uint64, error) {
	return nil, nil
}

// Returns h, v_i.
func (*Prover) createMaskedFunction() (*ring.Poly, []*ring.Poly, error) {
	return nil, nil, nil
}

// Returns c.
func (*Verifier) createChallenge() (*ring.Poly, error) {
	return nil, nil
}

// Returns z_i.
func (*Prover) createMaskedOpenings() ([][]*ring.Poly, error) {
	return nil, nil
}

// Completes the verification.
func (*Verifier) verify() (bool, error) {
	return false, nil
}
