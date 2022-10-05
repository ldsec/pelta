package main

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
	//"math"
	//"math/big"
	//"math/bits"
)

// Implementation of RLWE relations as a linear operation

func mainfunc() {

	// BFV parameters (128 bit security) with plaintext modulus 65929217
	paramDef := bfv.PN13QP218
	paramDef.T = 0x3ee0001

	params, err := bfv.NewParametersFromLiteral(paramDef)
	if err != nil {
		panic(err)
	}

	// Define the ring
	ringQP := params.RingQ()

	fmt.Printf("___ Implementation of RLWE relations ___\n\n")

	// DISCLAIMER
	/*
		TODO:
		-

		Note:
		-
	*/

	// Define PRNG
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	// Define Samplers
	// Uniform sampler
	uniformSampler := ring.NewUniformSampler(prng, params.RingQ())
	_ = uniformSampler

	// Ternary sampler  Montgomerry->False
	ternarySampler := ring.NewTernarySampler(prng, params.RingQ(), 1.0/3.0, false)
	_ = ternarySampler

	// Gaussian sampler
	gaussianSampler := ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), 16)
	_ = gaussianSampler

	// Define kgen
	kgen := bfv.NewKeyGenerator(params)

	// Define sk
	sk := kgen.GenSecretKey()
	s_poly := sk.Value.Q.CopyNew()
	s := sk.Value.Q.CopyNew()

	// Retrieve sk in Poly form
	ringQP.InvNTT(s_poly, s_poly)
	ringQP.InvMForm(s_poly, s_poly)

	// Define pk
	pk := kgen.GenPublicKey(sk)
	p0 := pk.Value[0].Q
	p1 := pk.Value[1].Q

	// [a] = pk[1]
	a := p1.CopyNew()
	_ = a

	// [e] = pk[0] + pk[1].s
	e := p0.CopyNew()
	ringQP.MulCoeffsMontgomeryAndAdd(a, s, e)

	tmp := e.CopyNew()
	ringQP.Neg(p1, a)
	ringQP.MulCoeffsMontgomeryAndAdd(a, s, tmp)

	if tmp.Equals(p0) == false {
		fmt.Print("tmp != p0 \n")
	} else {
		fmt.Printf("Reconstruction of p0 ok\n ")
	}

	e_poly := p0.CopyNew()
	ringQP.InvNTT(e, e_poly)
	ringQP.InvMForm(e_poly, e_poly)
	//fmt.Printf("e : %v\n", e_poly.Coeffs[0][0:20])

	// Reconstruction
	s_ntt := ringQP.NewPoly()
	ringQP.MForm(s_poly, s_ntt)
	ringQP.NTT(s_ntt, s_ntt)

	e_ntt := ringQP.NewPoly()
	ringQP.MForm(e_poly, e_ntt)
	ringQP.NTT(e_ntt, e_ntt)

	// p0 = -p1°s + e
	new_p0 := e_ntt.CopyNew()
	ringQP.MulCoeffsMontgomeryAndAdd(a, s_ntt, new_p0)

	if new_p0.Equals(p0) == false {
		fmt.Print("[FAIL]new_p0 != p0 in NTT and MForm \n")
	} else {
		fmt.Print("[OK] new_p0 != p0 in NTT and MForm \n")
	}

	ringQP.MForm(new_p0, new_p0)
	if new_p0.Equals(p0) == false {
		fmt.Print("[FAIL] MForm(new_p0) != p0 in NTT and MForm \n")
	} else {
		fmt.Print("[OK] MForm(new_p0) != p0 in NTT and MForm \n")
	}

}

func main() {
	mainfunc()
}
