package main

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"math"
	"math/big"
)

// Implementation of ENS20 automorphisms

func mainfunc() {

	// BFV parameters (128 bit security) with plaintext modulus 65929217
	paramDef := bfv.PN13QP218
	paramDef.T = 0x3ee0001

	params, err := bfv.NewParametersFromLiteral(paramDef)
	if err != nil {
		panic(err)
	}

	// Define the ring
	ringQP := params.RingQP()

	fmt.Printf("___ Implementation of ENS20 automorphism ___\n\n")

	// DISCLAIMER
	/*
		- TODO 1/k for all the different levels

		Note:
		-
	*/

	// Define PRNG
	prng, err := utils.NewPRNG()

	// Define Samplers
	uniformSampler := ring.NewUniformSampler(prng, params.RingQP())

	// Inputs
	m := ringQP.NewPoly()
	uniformSampler.Read(m)
	m.Coeffs[0][0] = 1

	// Define the parameters
	N := params.N()
	k := 4

	// Compute bigQ
	ModulusAtLevel := make([]*big.Int, len(ringQP.Modulus))
	ModulusAtLevel[0] = NewUint(ringQP.Modulus[0])
	for i := 1; i < len(ringQP.Modulus); i++ {
		val := new(big.Int).SetUint64(ringQP.Modulus[i])
		ModulusAtLevel[i] = new(big.Int).Mul(ModulusAtLevel[i-1], val)
	}
	//fmt.Printf("Q %v\n", ModulusAtLevel)

	// Invert k in Rq
	scalarBig := new(big.Int).ModInverse(new(big.Int).SetUint64(uint64(k)), ModulusAtLevel[len(ringQP.Modulus)-1])
	fmt.Printf("1/k:%v\n", scalarBig)

	// Check 1/k * k = 1 mod Q
	checkScalar := new(big.Int).Mul(scalarBig, new(big.Int).SetUint64(uint64(k)))
	checkScalarMod := new(big.Int).Mod(checkScalar, ModulusAtLevel[len(ringQP.Modulus)-1])
	fmt.Printf("1/k * k mod Q :%v\n", checkScalarMod)

	// Define the Galois Element 2N/k
	galEl := uint64(2*N/k + 1) //uint64(2*N/k+1)
	fmt.Printf("N:%v, k:%v, galEl:%v\n", N, k, galEl)
	//fmt.Printf("input in Poly -- %v\n", m.Coeffs[0][0:10])

	// Compute the automorphism sig for galEl:=2N/k+1
	mout := ringQP.NewPoly()
	ringQP.Permute(m, galEl, mout)
	//fmt.Printf("output in Poly-- %v\n", mout.Coeffs[0][0:10])

	// Compute the trace
	traceout := ringQP.NewPoly()
	for i := 0; i < k; i++ {
		ringQP.Permute(m, uint64(math.Pow(float64(galEl), float64(i))), mout)
		ringQP.Add(traceout, mout, traceout)
	}
	//fmt.Printf("traceout in Poly-- %v\n\n", traceout.Coeffs[0][0:10])

	// Test for multiple inputs
	m2 := ringQP.NewPoly()
	uniformSampler.Read(m2)
	m2.Coeffs[0][0] = 2

	m3 := ringQP.NewPoly()
	uniformSampler.Read(m3)
	m3.Coeffs[0][0] = 3

	m4 := ringQP.NewPoly()
	uniformSampler.Read(m4)
	m4.Coeffs[0][0] = 4

	/*
		fmt.Printf("m1 -- %v\n", m.Coeffs[0][0:10])
		fmt.Printf("m2 -- %v\n", m2.Coeffs[0][0:10])
		fmt.Printf("m3 -- %v\n", m3.Coeffs[0][0:10])
		fmt.Printf("m4 -- %v\n", m4.Coeffs[0][0:10])
	*/

	// Compute the trace
	traceout1 := ringQP.NewPoly()
	mout.Zero() //Cleaning the temporary value
	for i := 0; i < k; i++ {
		ringQP.Permute(m, uint64(math.Pow(float64(galEl), float64(i))), mout) // Executing sigma
		ringQP.Add(traceout1, mout, traceout1)                                // Adding to the trace
	}
	//fmt.Printf("traceout1 in Poly-- %v\n\n", traceout1.Coeffs[0][0:10])

	traceout2 := ringQP.NewPoly()
	mout.Zero()
	for i := 0; i < k; i++ {
		ringQP.Permute(m2, uint64(math.Pow(float64(galEl), float64(i))), mout)
		ringQP.Add(traceout2, mout, traceout2)
	}
	//fmt.Printf("traceout2 in Poly-- %v\n\n", traceout2.Coeffs[0][0:10])

	traceout3 := ringQP.NewPoly()
	mout.Zero()
	for i := 0; i < k; i++ {
		ringQP.Permute(m3, uint64(math.Pow(float64(galEl), float64(i))), mout)
		ringQP.Add(traceout3, mout, traceout3)
	}
	//fmt.Printf("traceout3 in Poly-- %v\n\n", traceout3.Coeffs[0][0:10])

	traceout4 := ringQP.NewPoly()
	mout.Zero()
	for i := 0; i < k; i++ {
		ringQP.Permute(m4, uint64(math.Pow(float64(galEl), float64(i))), mout)
		ringQP.Add(traceout4, mout, traceout4)
	}
	//fmt.Printf("traceout4 in Poly-- %v\n\n", traceout4.Coeffs[0][0:10])

	// Shift and sum the traces for i:=0..k-1
	traceout.Zero()
	ringQP.Shift(traceout1, 0, traceout1)
	ringQP.MulScalarBigint(traceout1, scalarBig, traceout1)
	ringQP.Add(traceout1, traceout, traceout)

	ringQP.Shift(traceout2, int(N-1), traceout2)
	ringQP.MulScalarBigint(traceout2, scalarBig, traceout2)
	ringQP.Add(traceout2, traceout, traceout)

	ringQP.Shift(traceout3, int(N-2), traceout3)
	ringQP.MulScalarBigint(traceout3, scalarBig, traceout3)
	ringQP.Add(traceout3, traceout, traceout)

	ringQP.Shift(traceout4, int(N-3), traceout4)
	ringQP.MulScalarBigint(traceout4, scalarBig, traceout4)
	ringQP.Add(traceout4, traceout, traceout)

	fmt.Printf("Lu k=4 -- %v\n\n", traceout.Coeffs[0][0:10])
	fmt.Printf("[WARNING] No 1/k yet\n")

}

func NewUint(v uint64) *big.Int {
	return new(big.Int).SetUint64(v)
}

func main() {
	mainfunc()
}
