package main

import (
	"testing"
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

func BenchmarkMain(b *testing.B){

	P := 10

	// BFV parameters (128 bit security) with plaintext modulus 65929217
	paramDef := bfv.PN13QP218
	paramDef.T = 0x3ee0001
	params, _ := bfv.NewParametersFromLiteral(paramDef)

	// Define the ring
	ringQ := params.RingQ()

	// Define list of polys per party
	polys := make([]*ring.Poly, P)

	// Enable random polynomial sampling
	prng, _ := utils.NewPRNG()
	uniformSamplerQ := ring.NewUniformSampler(prng, ringQ)

	// Generate a slice of polynomials polys[i] corresponds to party i's share
	for i:=0; i<P; i++{
		temp_poly := uniformSamplerQ.ReadNew()  

		for i := 0; i < params.N(); i++ {
			temp_poly.Coeffs[0][i] %= params.T()
		}
		ringQ.Reduce(temp_poly, temp_poly)
		polys[i] = temp_poly.CopyNew()
	}

	// Sum of the polynomials AggPoly(X) = sum_i P(X)
	AggPoly := ringQ.NewPoly()
	b.Run("sumPoly", func(b *testing.B) {
		for run := 0; run < b.N; run++ {

        	sumPoly(polys, AggPoly, ringQ, P)
        }
    })
	

	// Create the polynomial evaluations P(a) = sum_i P_i(a)
	point := uint64(2)

	evals := make([]uint64, P+1)  // Slice of the evaluations [P_1(a), ..., P_N(a), 0]
	for i:=0; i<P; i++{
		evals[i] = ring.EvalPolyModP(point, polys[i].Coeffs[0], params.T()) % params.T()
	}


	// Compute the aggregation on the scalars
	res := uint64(0)
	b.Run("Agg", func(b *testing.B) {
		for run := 0; run < b.N; run++ {
			res = sumEvals(res, AggPoly, evals, point, P, params)
		}
	})


}
