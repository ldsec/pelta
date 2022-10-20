package main

import (
	"fmt"

	"math"

	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
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

	//fmt.Printf("___ Implementation of Ajtai relations ___\n\n")
	q := params.Q()[0]
	fmt.Printf("q %v\n", q)
	fmt.Printf("     logq %v\n", int64(math.Log2(float64(q))))

	// Define PRNG
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	// Define Samplers
	// >> Uniform sampler
	uniformSampler := ring.NewUniformSampler(prng, params.RingQ())
	_ = uniformSampler

	// >> Ternary sampler  Montgomerry->False
	ternarySampler := ring.NewTernarySampler(prng, params.RingQ(), 1.0/3.0, false)
	_ = ternarySampler

	// Ajtai Commitment Scheme: c = A1.m + A2.r mod p
	// >> define p
	p := uint64(1<<20 - 9) // (1<<16 - 15)
	fmt.Printf("p %v\n", p)
	fmt.Printf("     logp %v\n\n", (math.Log2(float64(p))))

	// Get A1  of size 1 x N
	A1 := ringQP.NewPoly()
	uniformSampler.Read(A1)
	for i := 0; i < params.N(); i++ {
		for j, _ := range params.Q() {
			A1.Coeffs[j][i] = A1.Coeffs[0][i] % p
		}
	}
	//fmt.Printf("A1 %v\n", A1.Coeffs[0][0:10])

	// Get A2 of size 1 x N
	A2 := ringQP.NewPoly()
	uniformSampler.Read(A2)
	for i := 0; i < params.N(); i++ {
		for j, _ := range params.Q() {
			A2.Coeffs[j][i] = A2.Coeffs[0][i] % p
		}
	}
	//fmt.Printf("A2 %v\n", A2.Coeffs[0][0:10])

	fmt.Printf("[WARNING] A1 and A2 are over 0, p-1\n")

	// Sample message m
	m := ringQP.NewPoly()
	m_neg := ringQP.NewPoly()
	m_pos := ringQP.NewPoly()
	ternarySampler.Read(m)
	//m.Zero()
	for i := 0; i < params.N(); i++ {
		if m.Coeffs[0][i] > 3 {
			m_neg.Coeffs[0][i] = 1
			m_pos.Coeffs[0][i] = 0
		} else {
			m_pos.Coeffs[0][i] = m.Coeffs[0][i]
		}
	}
	//fmt.Printf("m  %v\n", m.Coeffs[0][0:10])

	// Sample randomness r
	r := ringQP.NewPoly()
	r_neg := ringQP.NewPoly()
	r_pos := ringQP.NewPoly()
	ternarySampler.Read(r)
	//r.Zero()
	for i := 0; i < params.N(); i++ {
		if r.Coeffs[0][i] > 3 {
			r_neg.Coeffs[0][i] = 1
			r_pos.Coeffs[0][i] = 0
		} else {
			r_pos.Coeffs[0][i] = r.Coeffs[0][i]

		}
	}

	//fmt.Printf("r  %v\n", r.Coeffs[0][0:10])

	///////////////////////////////////////////////////////////////////////////////////////////////
	fmt.Printf("\n-- Commitment with int64 (-1, 0, 1) level 0 --\n")
	// Create the commitment in Z by replacing the -1 in the m,r
	com := int64(0)

	for i := 0; i < params.N(); i++ {
		if m.Coeffs[0][i] == 1 {
			com = (com + int64(A1.Coeffs[0][i])*int64(m.Coeffs[0][i])) % int64(q)
		} else if m.Coeffs[0][i] > 3 {
			com = (com + int64(A1.Coeffs[0][i])*int64(-1)) % int64(q)
		}

		if r.Coeffs[0][i] == 1 {
			com = (com + int64(A2.Coeffs[0][i])*int64(r.Coeffs[0][i])) % int64(q)
		} else if r.Coeffs[0][i] > 3 {
			com = (com + int64(A2.Coeffs[0][i])*int64(-1)) % int64(q)
		}
	}

	k := int64(math.Round(float64(com-(com%int64(p))) / float64(p)))

	comp := com % int64(p)

	fmt.Printf(" [com]   level 0  %v\n", com)
	//fmt.Printf(" [com]q  level 0  %v\n", com%int64(q))
	//fmt.Printf(" [q+com]q level 0  %v\n", (int64(q)+com)%int64(q)) DON'T DO IT -- WRAPAROUND
	fmt.Printf(" [com]p  level 0  %v\n", comp)

	if comp != com-k*int64(p) {
		fmt.Printf("		k  %v\n", k)
		fmt.Printf("		[com] - k*p %v\n", com-k*int64(p))
		panic("[FAIL] mismatch [com]p != [com] - k*p \n")
	}

	if com != com%int64(q) {
		panic("[FAIL] mismatch [comp]q != [com] :: wraparound mod q\n")
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//fmt.Printf("\n-- Commitment with int64 (-1, 0, 1) level 1 --\n")
	// Create the commitment for level 1
	com3 := int64(0)

	for i := 0; i < params.N(); i++ {
		if m.Coeffs[0][i] == 1 {
			com3 = (com3 + int64(A1.Coeffs[0][i])*int64(m.Coeffs[0][i])) % int64(params.Q()[1])
		} else if m.Coeffs[0][i] > 3 {
			com3 = (com3 + int64(A1.Coeffs[0][i])*int64(-1)) % int64(params.Q()[1])
		}

		if r.Coeffs[0][i] == 1 {
			com3 = (com3 + int64(A2.Coeffs[0][i])*int64(r.Coeffs[0][i])) % int64(params.Q()[1])
		} else if r.Coeffs[0][i] > 3 {
			com3 = (com3 + int64(A2.Coeffs[0][i])*int64(-1)) % int64(params.Q()[1])
		}
	}
	//fmt.Printf(" [com]   level 1  %v\n", com3)
	//fmt.Printf(" [com]q  level 1  %v\n", com3%int64(params.Q()[1]))
	//fmt.Printf(" [com]p  level 1  %v\n", com3%int64(p))
	if com3%int64(p) != com%int64(p) {
		panic("[FAIL] mismatch [comp]p level 0 != [com]p level 1\n")
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//fmt.Printf("\n-- Commitment with int64 (0, 0, 1)-(1, 0, 0) level 0 --\n")
	// Create the commitment for level 1
	com4 := int64(0)

	tmp_res := ringQP.NewPoly()
	ringQP.MulCoeffs(A2, r_pos, tmp_res)
	ringQP.MulCoeffsAndAdd(A1, m_pos, tmp_res)

	tmp_res_neg := ringQP.NewPoly()
	ringQP.MulCoeffs(A2, r_neg, tmp_res_neg)
	ringQP.MulCoeffsAndAdd(A1, m_neg, tmp_res_neg)

	ringQP.Sub(tmp_res, tmp_res_neg, tmp_res)
	//fmt.Printf("++ %v\n", tmp_res.Coeffs[0][0:10])

	for i := 0; i < params.N(); i++ {
		com4 = (com4 + int64(tmp_res.Coeffs[0][i])%int64(params.Q()[0])) % int64(params.Q()[0])
	}

	//fmt.Printf(" [com]            %v\n", com4)
	//fmt.Printf(" [com]p           %v\n", com4%int64(p))

	///////////////////////////////////////////////////////////////////////////////////////////////
	//fmt.Printf("\n-- Commitment with uint64 (0, 0, 1)-(1, 0, 0) level 0 --\n")
	// Create the commitment for level 1
	com44 := uint64(0)

	for i := 0; i < params.N(); i++ {
		com44 = (com44 + uint64(tmp_res.Coeffs[0][i])%uint64(params.Q()[0])) % uint64(params.Q()[0])
	}

	//fmt.Printf(" [com]            %v\n", com44)
	//fmt.Printf(" [com]p           %v\n", com44%p)

	///////////////////////////////////////////////////////////////////////////////////////////////
	fmt.Printf("\n-- Commitment with uint64 (q-1, 0, 1) level 0 --\n")
	// Create the commitment in Z by replacing the -1 in the m,r
	com2 := uint64(0)

	for i := 0; i < params.N(); i++ {
		com2 = (com2 + uint64(A1.Coeffs[0][i])*uint64(m.Coeffs[0][i])) % uint64(params.Q()[0])
		com2 = (com2 + uint64(A2.Coeffs[0][i])*uint64(r.Coeffs[0][i])) % uint64(params.Q()[0])
	}
	fmt.Printf(" [com]    level 0  %v\n", com2)
	fmt.Printf(" [com]q   level 0  %v\n", (com2)%uint64(q))

	///////////////////////////////////////////////////////////////////////////////////////////////
	//fmt.Printf("\n-- Commitment with int64 (q-1, 0, 1) level 0 --\n")
	// Create the commitment for level 1
	com5 := int64(0)

	for i := 0; i < params.N(); i++ {
		com5 = (com5 + int64(A1.Coeffs[0][i])*int64(m.Coeffs[0][i])) % int64(params.Q()[0])
		com5 = (com5 + int64(A2.Coeffs[0][i])*int64(r.Coeffs[0][i])) % int64(params.Q()[0])
	}
	//fmt.Printf(" [com]   level 1  %v\n", com5)
	//fmt.Printf(" [com]q  level 1  %v\n", com5%int64(params.Q()[0]))
	//fmt.Printf(" [com]p  level 1  %v\n", com5%int64(p))

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	// Test with (0, 1, 2)
	for i := 0; i < params.N(); i++ {
		m.Coeffs[0][i] += 1
	}
	ringQP.Reduce(m, m)

	// Test with (0, 1, 2)
	for i := 0; i < params.N(); i++ {
		r.Coeffs[0][i] += 1
	}
	ringQP.Reduce(r, r)

	///////////////////////////////////////////////////////////////////////////////////////////////
	//fmt.Printf("\n-- Commitment with uint64 (0, 1, 2) level 0 --\n")
	// Create the commitment for level 1
	com012 := uint64(0)

	for i := 0; i < params.N(); i++ {
		com012 = (com012 + uint64(A1.Coeffs[0][i])*uint64(m.Coeffs[0][i])) % uint64(params.Q()[0])

		com012 = (com012 + uint64(A2.Coeffs[0][i])*uint64(r.Coeffs[0][i])) % uint64(params.Q()[0])
	}

	//fmt.Printf(" [com]   level 1  %v\n", com012)
	//fmt.Printf(" [com]q  level 1  %v\n", com012%(params.Q()[1]))
	//fmt.Printf(" [com]p  level 1  %v\n", com012%(p))

	///////////////////////////////////////////////////////////////////////////////////////////////
	//fmt.Printf("\n-- Commitment with uint64 (0, 1, 2) level 1 --\n")
	// Create the commitment for level 1
	com012_ := uint64(0)

	for i := 0; i < params.N(); i++ {
		com012_ = (com012_ + uint64(A1.Coeffs[0][i])*uint64(m.Coeffs[0][i])) % uint64(params.Q()[1])

		com012_ = (com012_ + uint64(A2.Coeffs[0][i])*uint64(r.Coeffs[0][i])) % uint64(params.Q()[1])
	}

	//fmt.Printf(" [com]   level 1  %v\n", com012_)
	//fmt.Printf(" [com]q  level 1  %v\n", com012_%(params.Q()[1]))
	//fmt.Printf(" [com]p  level 1  %v\n", com012_%(p))
	if com012%(params.Q()[1]) != com012_%(params.Q()[1]) {
		panic("[FAIL] mismatch [com]p l0 != [com]p l1 for (0,1,2)\n")
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	// CAN WE RECONSTRUCT m AND r?

	ones := ringQP.NewPoly()
	for i := 0; i < params.N(); i++ {
		ones.Coeffs[0][i] = p - 1 //18014398508400640
	}

	ringQP.Add(m, ones, m)
	ringQP.Add(r, ones, r)
	//fmt.Printf("m  %v\n", m.Coeffs[0][0:10])
	//fmt.Printf("r  %v\n", r.Coeffs[0][0:10])

	///////////////////////////////////////////////////////////////////////////////////////////////
	fmt.Printf("\n-- Commitment with uint64 (0, 1, 2)-(1, 1, 1) level 0 --\n")
	// Create the commitment in Z by replacing the -1 in the m,r
	com_rec := uint64(0)

	for i := 0; i < params.N(); i++ {
		com_rec = (com_rec + uint64(A1.Coeffs[0][i])*uint64(ones.Coeffs[0][i])) % uint64(params.Q()[0])
		com_rec = (com_rec + uint64(A2.Coeffs[0][i])*uint64(ones.Coeffs[0][i])) % uint64(params.Q()[0])
	}
	//fmt.Printf(" [com]    level 0  %v\n", com_rec)
	//fmt.Printf(" [com]q   level 0  %v\n", (com_rec)%uint64(q))
	com_rec = com_rec + com012
	fmt.Printf(" [com]    level 0  %v\n", (com_rec))
	//fmt.Printf(" [com]q   level 0  %v\n", (com_rec)%uint64(q))
	fmt.Printf(" [com]p   level 0  %v\n", (com_rec)%uint64(p))

	if (com_rec) != (com_rec)%uint64(q) {
		panic("[FAIL] mismatch [com] != [com]q :: wraparound mod q\n")
	}

	if ((p + com_rec) % uint64(p)) != uint64((int64(p)+comp)%int64(p)) {
		fmt.Printf("[com -1,0,1]p       %v\n", uint64((int64(p)+comp)%int64(p)))
		fmt.Printf("[com 0,1,2]p - f(1) %v\n", ((p + com_rec) % uint64(p)))
		panic("[FAIL] mismatch [com -1,0,1]p != [com 0,1,2]p - (A1.1+A2.2)\n")
	}

	k2 := uint64(math.Round(float64(com_rec-((com_rec)%p)) / float64(p)))

	if (com_rec)%p != com_rec-k2*p {
		fmt.Printf("		k  %v\n", k2)
		fmt.Printf("		[com] - k*p %v\n", com_rec-k2*p)
		panic("[FAIL] mismatch [com]p != [com] - k*p \n")
	} else {
		fmt.Printf("Success\n")
	}

}

func main() {
	mainfunc()
}
