package main

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/bfv"
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

	// Define kgen
	kgen := bfv.NewKeyGenerator(params)

	// Define sk
	sk := kgen.GenSecretKey()
	s_poly := sk.Value.Q.CopyNew()
	s := sk.Value.Q.CopyNew()

	// Retrieve sk in Poly form
	ringQP.InvNTT(s_poly, s_poly)
	ringQP.InvMForm(s_poly, s_poly)
	fmt.Printf("s : %v\n", s_poly.Coeffs[0][0:20])

	// Define pk
	pk := kgen.GenPublicKey(sk)
	p0 := pk.Value[0].Q
	p1 := pk.Value[1].Q

	// [a] = pk[1]
	p1_poly := p1.CopyNew()
	ringQP.InvNTT(p1_poly, p1_poly)
	ringQP.InvMForm(p1_poly, p1_poly)

	// [e] = pk[0] + pk[1].s
	e := p0.CopyNew()
	ringQP.MulCoeffsMontgomeryAndAdd(p1, s, e)

	e_poly := p0.CopyNew()
	ringQP.InvNTT(e, e_poly)
	ringQP.InvMForm(e_poly, e_poly)
	fmt.Printf("e : %v\n", e_poly.Coeffs[0][0:20])

	// Reconstruction
	s_ntt := ringQP.NewPoly()
	ringQP.NTT(s_poly, s_ntt)
	ringQP.MForm(s_ntt, s_ntt)

	e_ntt := ringQP.NewPoly()
	ringQP.NTT(e_poly, e_ntt)
	ringQP.MForm(e_ntt, e_ntt)

	// p0 = -p1°s + e
	neg_a := p1.CopyNew()
	ringQP.Neg(p1, neg_a)

	new_p0 := e_ntt.CopyNew()
	ringQP.MulCoeffsMontgomeryAndAdd(neg_a, s_ntt, new_p0)

	if new_p0.Equals(p0) == false {
		fmt.Print("[FAIL]new_p0 != p0 in NTT and MForm \n")
	} else {
		fmt.Print("[OK] new_p0 != p0 in NTT and MForm \n")
	}

	// Reconstruction without MForm
	s_ntt_noM := ringQP.NewPoly()
	ringQP.NTT(s_poly, s_ntt_noM)

	e_ntt_noM := ringQP.NewPoly()
	ringQP.NTT(e_poly, e_ntt_noM)

	a_ntt_noM := ringQP.NewPoly()
	ringQP.NTT(p1_poly, a_ntt_noM)

	ringQP.Neg(a_ntt_noM, a_ntt_noM)

	new_p0_noM := e_ntt_noM.CopyNew()
	ringQP.MulCoeffsAndAdd(a_ntt_noM, s_ntt_noM, new_p0_noM)

	ringQP.MForm(new_p0_noM, new_p0_noM)

	if new_p0_noM.Equals(p0) == false {
		fmt.Print("[FAIL] MForm(new_p0) != p0 in NTT and MForm \n")
	} else {
		fmt.Print("[OK] MForm(new_p0) != p0 in NTT and MForm \n")
	}

	/////////////////////////////////////////////////////////////////////////
	// Reconstruction with sk in (0,1,2)
	s_poly2 := s_poly.CopyNew()
	for i := 0; i < params.N(); i++ {
		s_poly2.Coeffs[0][i] += 1
	}
	ringQP.Reduce(s_poly2, s_poly2)
	fmt.Printf("s : %v\n", s_poly2.Coeffs[0][0:20])

	ones := ringQP.NewPoly()
	for i := 0; i < params.N(); i++ {
		ones.Coeffs[0][i] = 1
	}

	// Encyption of ones without MForm
	ringQP.NTT(ones, ones)
	onesEnc := ringQP.NewPoly()
	ringQP.MulCoeffs(a_ntt_noM, ones, onesEnc)

	// Encryption of sk in (012)
	s_ntt_012 := ringQP.NewPoly()
	ringQP.NTT(s_poly2, s_ntt_012)

	p0_012 := e_ntt_noM.CopyNew()
	ringQP.MulCoeffsAndAdd(a_ntt_noM, s_ntt_012, p0_012)

	// Sub p0_012-oneEnc
	ringQP.Sub(p0_012, onesEnc, p0_012)

	// MForm
	ringQP.MForm(p0_012, p0_012)
	if p0_012.Equals(p0) == false {
		fmt.Print("[FAIL] MForm(p0_012) != p0 in NTT and MForm \n")
	} else {
		fmt.Print("[OK] MForm(p0_012) != p0 in NTT and MForm \n")
	}

}

func main() {
	mainfunc()
}
