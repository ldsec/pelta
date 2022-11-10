package fastens20

import (
	"github.com/ldsec/codeBase/commitment/crypto"
	"github.com/ldsec/codeBase/commitment/fastmath"
)

// PublicParams contains the public parameters of the protocol.
type PublicParams struct {
	config crypto.Config
	A      *fastmath.IntMatrix
	U      *fastmath.IntVec
	B0     *fastmath.PolyNTTMatrix
	B      *fastmath.PolyNTTMatrix
	Sig    fastmath.Automorphism
}

func NewPublicParameters(sis crypto.SISProblem, config crypto.Config) PublicParams {
	bSize := config.NumSplits() + 3
	B0 := fastmath.NewRandomPolyMatrix(config.Kappa, config.Lambda+config.Kappa+bSize, config.UniformSampler, config.BaseRing)
	b := fastmath.NewRandomPolyMatrix(bSize, B0.Cols(), config.UniformSampler, config.BaseRing)
	sig := fastmath.NewAutomorphism(uint64(config.D), uint64(config.K))
	return PublicParams{
		config: config,
		A:      sis.A,
		U:      sis.U,
		B0:     B0.NTT(),
		B:      b.NTT(),
		Sig:    sig,
	}
}
