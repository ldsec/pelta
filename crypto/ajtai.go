package crypto

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/tuneinsight/lattigo/v4/ring"
)

type AjtaiCommitment struct {
	A    *fastmath.IntMatrix
	B    *fastmath.IntMatrix
	s    *fastmath.IntVec
	r    *fastmath.IntVec
	comQ *fastmath.IntVec
	p    *big.Int
}

func NewAjtaiCommitment(s, r *fastmath.IntVec, comSize int, p *big.Int, q *big.Int, baseRing *ring.Ring) AjtaiCommitment {
	A := fastmath.NewRandomIntMatrix(s.Size(), comSize, p, baseRing)
	B := fastmath.NewRandomIntMatrix(s.Size(), comSize, p, baseRing)
	// Compute kappa.
	comQ := A.MulVec(s).Add(B.MulVec(r))
	comPNeg := comQ.Copy().Reduce(p).Neg()
	// comQ - comP = (As + Br) - [(As + Br) mod p]
	diff := comQ.Copy().Add(comPNeg)
	// kappa := (comQ - comP) / p s.t. kappa * p is the difference between
	kappa := fastmath.NewIntVec(diff.Size(), baseRing)
	kappa.Populate(func(i int) uint64 {
		return diff.Get(i) / p.Uint64()
	})
	return AjtaiCommitment{A, B, s, r, comQ, p}
}