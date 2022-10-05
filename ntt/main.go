package main

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/ldsec/lattigo/v2/bfv"
)

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

	// Define PRNG
	prng, err := utils.NewPRNG()

	// Ternary sampler  Montgomerry->False
	ternarySampler := ring.NewTernarySampler(prng, params.RingQP(), 1.0/3.0, false)
	_ = ternarySampler
	
	// Uniform sampler
	uniformSampler := ring.NewUniformSampler(prng, params.RingQP())
	_ = uniformSampler



// Creating a random polynomial
	m:=ringQP.NewPoly()
	uniformSampler.Read(m)
	//m.Zero()
	//m.Coeffs[0][0]=1
	fmt.Printf("     m : %v\n", m.Coeffs[0][0:10])
	//      m : [1 0 0 0 0 0 0 0 0 0]

// Computing its NTT transform
	refPoly := ringQP.NewPoly()
	ringQP.NTT(m, refPoly)
	fmt.Printf("NTT(m) : %v\n", refPoly.Coeffs[0][0:10])
	//NTT(m) : [1 1 1 1 1 1 1 1 1 1]


// Constructing the NTT matrix !!!!!! NOTE: only consider the first level q0 !!!!!!

	w := ring.InvMForm(ringQP.NttPsi[0][ringQP.N>>1], ringQP.Modulus[0], ringQP.MredParams[0])


	//T[i, j] = w^{j*(2*rev(i)+1)} mod q  represented as N poly corresponding to each row of T.
	T := make([]*ring.Poly, ringQP.N)
	mask := uint64(2*ringQP.N-1)

	for i := 0; i < ringQP.N; i++ {
		T[i] = ringQP.NewPoly()

		twoirev := 2*utils.BitReverse64(uint64(i), uint64(params.LogN()))+1

		for j := 0; j < ringQP.N; j++{

			gen := uint64(j) * twoirev & mask

			T[i].Coeffs[0][j] = ring.ModExp(w, int(gen), ringQP.Modulus[0])
		}
	}


	// Create T.vec(s)
	u := ringQP.NewPoly()

	for i:=0; i<ringQP.N; i++{
		// get a temporary vector to store T[i]°s
		tmp := ringQP.NewPoly()

		// Mult s in poly form by T[i] coeff-wise
		ringQP.MulCoeffs(m, T[i], tmp)

		// Inner prod to get \sum_{j=0}^{N-1} w^{ij} * s_{j} in the 1-th slot of tmp
		for j := 0; j < params.LogN(); j++ {
			tmp2 := ringQP.NewPoly()
			ringQP.Shift(tmp, 1<<j, tmp2)
			ringQP.Add(tmp, tmp2, tmp)
		}

		// Copy the 1-st slot of tmp into the i-th slot of u.
		u.Coeffs[0][i] = tmp.Coeffs[0][0]
	}

fmt.Printf("  m x T: %10d\n", u.Coeffs[0][0:10])

cp_u := ringQP.NewPoly()
cp_ref := ringQP.NewPoly()
for i:=0; i<ringQP.N; i++{
	cp_u.Coeffs[0][i] = u.Coeffs[0][i]
	cp_ref.Coeffs[0][i] = refPoly.Coeffs[0][i]
} 

if cp_u.Equals(cp_ref) == false {
	panic("[FAIL] NTT not verified")
}else{
	fmt.Printf("[PASS] NTT matrix verified\n")
}


}


func main() {
	mainfunc()
}
