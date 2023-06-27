package main

import (
	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
	"log"
	"time"
	"fmt"
)

func sumPoly(polys []*ring.Poly, AggPoly *ring.Poly, ringQ *ring.Ring, P int) {
	for i:=0; i<P; i++{
		ringQ.Add(AggPoly, polys[i], AggPoly)
		ringQ.Reduce(AggPoly, AggPoly)
	}
}


func sumEvals(res uint64, AggPoly *ring.Poly, evals []uint64, point uint64, P int, params bfv.Parameters) uint64 {

	// Store AggPoly(a) at the end of the slice: i.e., [P_1(a), ..., P_N(a), AggPoly(a)]
	evals[P] = ring.EvalPolyModP(point, AggPoly.Coeffs[0], params.T()) % params.T()

	// Sum the scalars
	for i:=0; i<P; i++{
		res = (res + evals[i]) % params.T()
	}
	return res
}

func mainfunc(P int) (time.Duration, time.Duration){

	fmt.Printf("Aggregation of %v parties\n", P)

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
	start := time.Now()
	sumPoly(polys, AggPoly, ringQ, P)
	elapsed := time.Since(start)
	log.Printf("Poly Sum: %s", elapsed)

	// Create the polynomial evaluations P(a) = sum_i P_i(a)
	point := uint64(2)

	evals := make([]uint64, P+1)  // Slice of the evaluations [P_1(a), ..., P_N(a), 0]
	for i:=0; i<P; i++{
		evals[i] = ring.EvalPolyModP(point, polys[i].Coeffs[0], params.T()) % params.T()
	}



	// Compute the aggregation on the scalars
	res := uint64(0)
	startAgg := time.Now()
	res = sumEvals(res, AggPoly, evals, point, P, params)
	elapsedAgg := time.Since(startAgg)
	log.Printf("Agg  Sum: %s", elapsedAgg)


	if res != evals[P] {
		panic("Aggregation mismatch")
	}
	return elapsed, elapsedAgg
}

func main() {
	mainfunc(2)
}