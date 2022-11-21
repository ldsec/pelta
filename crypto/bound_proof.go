package crypto

import (
	"errors"
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
)

type BoundProofParams struct {
	B *fastmath.IntMatrix
	y *fastmath.IntVec
}

func NewBoundProofParams(rows, cols int, n *big.Int, baseRing *ring.Ring) BoundProofParams {
	B := fastmath.NewRandomIntMatrix(rows, cols, n, baseRing)
	y := fastmath.NewRandomIntVec(rows, n, baseRing)
	return BoundProofParams{B, y}
}

type BoundProof struct {
	BoundProofParams
	z *fastmath.IntVec
}

// NewBoundProof constructs a new bound proof
func NewBoundProof(x *fastmath.IntVec, rejSamplingBound uint64, params BoundProofParams) (BoundProof, error) {
	z := params.B.MulVec(x).Add(params.y)
	if !fastmath.AcceptIntSample(z, rejSamplingBound) {
		return BoundProof{}, errors.New("rejection sampling")
	}
	return BoundProof{params, z}, nil
}

// Verify verifies the bound proof against the given bound.
func (p BoundProof) Verify(bound uint64) bool {
	return p.z.Max() <= bound
}
