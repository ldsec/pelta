package main

import (
	"fmt"
	"github.com/tuneinsight/lattigo/v4/bfv"
	//"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Implementation of ENS20 commitments

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

	fmt.Printf("___ Implementation of ENS20 ___\n\n")

	// DISCLAIMER: Sep 23  SC
	/*
	   	- The sigma automorphism is not yet implemneted so we only work with k=1 for now
	    	- The shift 1/k*X^u is not implemented. It can be achieved with a Shift in lattigo
	    	- The inverse 1/k in Zq not implemented -- for now 1
	    	- The Abort and bound check in infinity norm are not implemented yet

	    	Note:
	    	- All values are in the Poly domain (except Atgammak)
	    	- Beaware of the changes in the protocol e.g., the NTT/iNTT operations -- They could be factored to speed up the code
	*/

	// Initialisation
	// SetElements the commitment params
	N := uint64(params.N()) // dimension of Rq
	n := 1                  // dimension for MSIS
	m1 := 1                 // dimension of the input in Rq
	m2 := 1                 // dimension of the input in Rq
	lbd := 1                // dimension for MLWE
	sizeB := n + lbd + 2    // dimension of bi
	size_t := n + m1 + 1    // dimension of the commitment
	size_bVec := size_t - n // dimension of bVec vector of {bi}
	k := 1                  // repeat rate in ENS20
	d1 := 16                // delta1 bound on the mask // Change for appropriate value

	// Define PRNG
	prng, err := utils.NewPRNG()

	// Define Samplers
	// Test uniform sampler
	uniformSampler := ring.NewUniformSampler(prng, params.RingQP())

	// Test ternary sampler  Montgomerry->False
	ternarySampler := ring.NewTernarySampler(prng, params.RingQP(), 1.0/3.0, false)

	// Gaussian sampler
	gaussianSampler := ring.NewGaussianSampler(prng, params.RingQP(), params.Sigma(), d1)

	// Inputs
	m := make([]*ring.Poly, m1)
	for i := 0; i < m1; i++ {
		m[i] = ringQP.NewPoly()
		uniformSampler.Read(m[i])
		// v DEBUG v -- SetElements the input vector
		//m[i].Zero()
	}

	// Create the corresponding polynomial
	mHat := make([]*ring.Poly, m1)
	for i := 0; i < m1; i++ {
		mHat[i] = ringQP.NewPoly()
		ringQP.InvNTT(m[i], mHat[i])
	}

	var A = make([][]*ring.Poly, m2)
	for i := 0; i < m2; i++ {
		A[i] = make([]*ring.Poly, m1)
		for j := 0; j < m1; j++ {
			A[i][j] = ringQP.NewPoly()
			uniformSampler.Read(A[i][j])
			// v DEBUG v -- SetElements the combination matrix A
			//A[i][j].Zero()
		}
	}

	var u = make([]*ring.Poly, m2)
	for i := 0; i < m2; i++ {
		u[i] = ringQP.NewPoly()
		for j := 0; j < m1; j++ {
			ringQP.MulCoeffsAndAdd(A[i][j], m[j], u[i])
		}
	}

	// Sample the public elements B in Rq^(n x (n+lbd+2)) in Poly
	var B = make([][]*ring.Poly, n)
	for i := 0; i < n; i++ {
		B[i] = make([]*ring.Poly, sizeB)
		for j := 0; j < sizeB; j++ {
			B[i][j] = ringQP.NewPoly()
			uniformSampler.Read(B[i][j])
		}
	}

	// Sample the public elements bVec in Rq^(n+lbd+2)xsizebVec in Poly
	var bVec = make([][]*ring.Poly, size_bVec)
	for i := 0; i < size_bVec; i++ {
		bVec[i] = make([]*ring.Poly, sizeB)
		for j := 0; j < sizeB; j++ {
			bVec[i][j] = ringQP.NewPoly()
			uniformSampler.Read(bVec[i][j])
		}
	}

	// Prover samples
	// Sample a polynomial g s.t. g_0=...=g_(k-1)=0
	g := ringQP.NewPoly()
	uniformSampler.Read(g)
	for i := 0; i < k; i++ {
		g.Coeffs[0][i] = 0
	}

	// Sample the commitment randomness r in Rq^((lbd+n+2))
	var r = make([]*ring.Poly, sizeB)
	for i := 0; i < sizeB; i++ {
		r[i] = ringQP.NewPoly()
		ternarySampler.Read(r[i])
		// v DEBUG v -- set randomness to zero
		//r[i].Zero()
	}

	// Create the commitments in Poly form
	var t = make([]*ring.Poly, size_t)

	// t0 = Br
	for i := 0; i < n; i++ {
		t[i] = ringQP.NewPoly()
		for j := 0; j < sizeB; j++ {
			ringQP.NTT(B[i][j], B[i][j])
			ringQP.NTT(r[j], r[j])
			ringQP.MulCoeffsAndAdd(B[i][j], r[j], t[i]) // Coefficient-wise multiplication
			ringQP.InvNTT(B[i][j], B[i][j])
			ringQP.InvNTT(r[j], r[j])
		}
		ringQP.InvNTT(t[i], t[i])
	}

	// t1[k] = <b1[k], r> + m^[k]  with m^ the invNTT transform of m
	for i := 0; i < m1; i++ {
		t[n+i] = ringQP.NewPoly()
		// create the inner product <b1, r>
		for j := 0; j < sizeB; j++ {
			ringQP.NTT(bVec[i][j], bVec[i][j])
			ringQP.NTT(r[j], r[j])
			ringQP.MulCoeffsAndAdd(bVec[i][j], r[j], t[n+i]) // coefficient-wise mult and accum in t1
			ringQP.InvNTT(bVec[i][j], bVec[i][j])
			ringQP.InvNTT(r[j], r[j])
		}
		ringQP.InvNTT(t[n+i], t[n+i])
		ringQP.Add(t[n+i], mHat[i], t[n+i]) //
	}

	// t2 = <b2, r> + g
	t[n+m1] = ringQP.NewPoly()
	// create the inner product <b1, r>  in Poly form
	for i := 0; i < sizeB; i++ {
		ringQP.NTT(bVec[m1][i], bVec[m1][i])
		ringQP.NTT(r[i], r[i])
		ringQP.MulCoeffsAndAdd(bVec[m1][i], r[i], t[n+m1]) // coefficient-wise mult and accum in t1
		ringQP.InvNTT(bVec[m1][i], bVec[m1][i])
		ringQP.InvNTT(r[i], r[i])
	}
	ringQP.InvNTT(t[n+m1], t[n+m1])
	ringQP.Add(t[n+m1], g, t[n+m1])

	// Create the mask
	var w = make([][]*ring.Poly, k)
	var y = make([][]*ring.Poly, k)
	for i := 0; i < k; i++ {
		w[i] = make([]*ring.Poly, n)
		y[i] = make([]*ring.Poly, sizeB)
		for j := 0; j < sizeB; j++ {
			y[i][j] = ringQP.NewPoly()
			gaussianSampler.Read(y[i][j])
			// v DEBUG v -- SetElements yi to zero
			//y[i][j].Zero()  // TODO rmv
		}

		for j := 0; j < n; j++ {
			w[i][j] = ringQP.NewPoly()
			for jj := 0; jj < sizeB; jj++ {
				ringQP.NTT(B[j][jj], B[j][jj])
				ringQP.NTT(y[i][jj], y[i][jj])
				ringQP.MulCoeffsAndAdd(B[j][jj], y[i][jj], w[i][j]) // Coefficient-wise multiplication
				ringQP.InvNTT(B[j][jj], B[j][jj])
				ringQP.InvNTT(y[i][jj], y[i][jj])
			}
			ringQP.InvNTT(w[i][j], w[i][j])
		}
	}

	// Generate challenge gamma_0...gamma_(k-1) Zq^m
	var gamma = make([][]*ring.Poly, k)
	for i := 0; i < k; i++ {
		gamma[i] = make([]*ring.Poly, m2)
		for j := 0; j < m2; j++ {
			gamma[i][j] = ringQP.NewPoly()
			uniformSampler.Read(gamma[i][j])
		}
	}

	// Create h
	// Create At
	At := transpose(A)
	temp := ringQP.NewPoly() // temp poly used through the script

	// Create the target value (At*gamma_k) -- in COEFF form
	Atgammak := make([][]*ring.Poly, k)
	for i := 0; i < k; i++ {
		Atgammak[i] = make([]*ring.Poly, m1)
		for j := 0; j < m1; j++ {
			Atgammak[i][j] = ringQP.NewPoly()
			for jj := 0; jj < m2; jj++ {
				ringQP.MulCoeffsAndAdd(At[j][jj], gamma[i][jj], Atgammak[i][j])
			}
		}
	}

	// Compute the inner product cste[i]=-<u, gamma[i]>  -- constant Poly
	cste := make([]*ring.Poly, k)
	for i := 0; i < k; i++ {
		cste[i] = ringQP.NewPoly()
		for j := 0; j < m2; j++ {
			ringQP.MulCoeffsAndAdd(u[j], gamma[i][j], cste[i])
		}

		for j := 0; j < params.LogN(); j++ {
			ringQP.Shift(cste[i], 1<<j, temp)
			ringQP.Add(cste[i], temp, cste[i])
		}
		tmpInner := cste[i].Coeffs[0][0]
		cste[i].Zero()
		cste[i].Coeffs[0][0] = tmpInner
		// Negate the inner product
		ringQP.Neg(cste[i], cste[i])
		ringQP.Reduce(cste[i], cste[i])
	}

	// Create  iNTT(N.Atgamma ° m[j]) = iNTT(N.Atgamma)* m^[j] -- in Poly form
	prodk := make([][]*ring.Poly, k)
	for i := 0; i < k; i++ {
		prodk[i] = make([]*ring.Poly, m1)
		for j := 0; j < m1; j++ {
			prodk[i][j] = ringQP.NewPoly()
			ringQP.MulCoeffs(Atgammak[i][j], m[j], prodk[i][j])
			ringQP.MulScalar(prodk[i][j], N, prodk[i][j])
			ringQP.InvNTT(prodk[i][j], prodk[i][j])
		}
	}

	// Sum_v=0^m1 [iNTT(N.Atgamma ° m_v)] - <u, gamma> -- In Poly form
	for i := 0; i < k; i++ {
		for j := 1; j < m1; j++ {
			ringQP.Add(prodk[i][0], prodk[i][j], prodk[i][0])
		}
		ringQP.Add(prodk[i][0], cste[i], prodk[i][0])
	}

	// Create the target value T = Sum_v=0^k-1 \sig_v[ prod_k ] -- In Poly form
	sumSig := make([]*ring.Poly, k)
	for i := 0; i < k; i++ {
		sumSig[i] = ringQP.NewPoly()
		// TODO sig is not a rotation!!!
		ringQP.Add(sumSig[i], prodk[i][0], sumSig[i])
	}

	// Create the value f = Sum_k X^k T[k]  -- In Poly form TODO
	fValue := ringQP.NewPoly()
	fValue = sumSig[0].CopyNew()
	for i := 0; i < k; i++ {
		// TODO
	}

	// h = g + f  -- In Poly form
	h := ringQP.NewPoly()
	ringQP.Add(g, fValue, h)

	// Create v_i
	// Create the value Atgamma_b1m[i][j][jj] = iNTT(N.Atgamma[i][j] ° b1[jj]) for jj in sizeB  -- in Poly form
	Atgamma_b1m := make([][][]*ring.Poly, k)
	for i := 0; i < k; i++ {
		Atgamma_b1m[i] = make([][]*ring.Poly, m1)
		for j := 0; j < m1; j++ {
			Atgamma_b1m[i][j] = make([]*ring.Poly, sizeB)
			for jj := 0; jj < sizeB; jj++ {
				Atgamma_b1m[i][j][jj] = ringQP.NewPoly()
				ringQP.NTT(bVec[j][jj], bVec[j][jj])
				ringQP.MulCoeffs(Atgammak[i][j], bVec[j][jj], Atgamma_b1m[i][j][jj])
				ringQP.MulScalar(Atgamma_b1m[i][j][jj], N, Atgamma_b1m[i][j][jj])
				ringQP.InvNTT(bVec[j][jj], bVec[j][jj])
				ringQP.InvNTT(Atgamma_b1m[i][j][jj], Atgamma_b1m[i][j][jj])
			}
		}
	}

	// Create the scalar product scalarProdVm = <Atgamma_b1m[0]
	scalarProdVm := make([][][][]*ring.Poly, k)
	for i := 0; i < k; i++ {
		scalarProdVm[i] = make([][][]*ring.Poly, k)
		for mu := 0; mu < k; mu++ {
			scalarProdVm[i][mu] = make([][]*ring.Poly, k)
			for nu := 0; nu < k; nu++ {
				scalarProdVm[i][mu][nu] = make([]*ring.Poly, m1)
				for j := 0; j < m1; j++ {
					scalarProdVm[i][mu][nu][j] = ringQP.NewPoly()
					for jj := 0; jj < sizeB; jj++ {
						ringQP.NTT(Atgamma_b1m[mu][j][jj], Atgamma_b1m[mu][j][jj])
						ringQP.NTT(y[(k+i-nu)%k][jj], y[(k+i-nu)%k][jj])
						ringQP.MulCoeffsAndAdd(Atgamma_b1m[mu][j][jj], y[(k+i-nu)%k][jj], scalarProdVm[i][mu][nu][j])
						ringQP.InvNTT(Atgamma_b1m[mu][j][jj], Atgamma_b1m[mu][j][jj])
						ringQP.InvNTT(y[(k+i-nu)%k][jj], y[(k+i-nu)%k][jj])
					}
					ringQP.InvNTT(scalarProdVm[i][mu][nu][j], scalarProdVm[i][mu][nu][j])
				}
			}
		}
	}

	// Sum the sum_mu sum_nu sum_j scalarProdVm[i][mu][nu][j]
	scalarProdV := make([]*ring.Poly, k)
	for i := 0; i < k; i++ {
		scalarProdV[i] = ringQP.NewPoly()

		for mu := 0; mu < k; mu++ {
			for nu := 0; nu < k; nu++ {
				for j := 0; j < m1; j++ {

					// TODO 1/k.X^mu

					// TODO sig^nu

					ringQP.Add(scalarProdV[i], scalarProdVm[i][mu][nu][j], scalarProdV[i])
				}
			}
		}
	}

	// Create the scalar product <b2, yi>
	scalarProd_b2y := make([]*ring.Poly, k)
	for i := 0; i < k; i++ {
		scalarProd_b2y[i] = ringQP.NewPoly()
		for j := 0; j < sizeB; j++ {
			ringQP.NTT(bVec[m1][j], bVec[m1][j])
			ringQP.NTT(y[i][j], y[i][j])
			ringQP.MulCoeffsAndAdd(bVec[m1][j], y[i][j], scalarProd_b2y[i])
			ringQP.InvNTT(bVec[m1][j], bVec[m1][j])
			ringQP.InvNTT(y[i][j], y[i][j])
		}
		ringQP.InvNTT(scalarProd_b2y[i], scalarProd_b2y[i])
	}

	// Create vi = scalarProdV + <b2, yi>
	vVector := make([]*ring.Poly, k)
	for i := 0; i < k; i++ {
		vVector[i] = ringQP.NewPoly()
		ringQP.Add(scalarProdV[i], scalarProd_b2y[i], vVector[i])
	}

	// Get challenge c in C: {-1, 0, 1}^N
	c := ringQP.NewPoly()
	ternarySampler.Read(c)
	// v DEBUG v -- SetElements challenge to vec(1) // TODO rmv
	//c.Zero()

	// Create the z_i openings
	var z = make([][]*ring.Poly, k)
	for i := 0; i < k; i++ {
		z[i] = make([]*ring.Poly, sizeB)

		for j := 0; j < sizeB; j++ {
			z[i][j] = ringQP.NewPoly()
			factor := ringQP.NewPoly()
			//ringQP.Shift // TODO sigma not a shift
			factor = c.CopyNew()
			ringQP.NTT(r[j], r[j])
			ringQP.NTT(factor, factor)
			ringQP.MulCoeffs(factor, r[j], z[i][j])
			ringQP.InvNTT(factor, factor)
			ringQP.InvNTT(r[j], r[j])
			ringQP.InvNTT(z[i][j], z[i][j])

			ringQP.Add(z[i][j], y[i][j], z[i][j])
		}
	}

	// Check the norm and the abort condition
	fmt.Printf("HOLD: TODO: add the abort condition on the norm\n")

	// Verify
	// Check the norm z_i
	fmt.Printf("HOLD: TODO: check norm zi\n")
	fmt.Printf("[CHECK] ||zi|| < Beta           -- NOT IMPLEMENTED TODO \n")

	// B.z =?= w + sigma(c)t0
	var t0_check = make([][]*ring.Poly, k)
	for i := 0; i < k; i++ {
		t0_check[i] = make([]*ring.Poly, n)

		tmpC := ringQP.NewPoly()

		// TODO sigma^i(c)
		tmpC = c.CopyNew()

		ringQP.NTT(tmpC, tmpC)

		for j := 0; j < n; j++ {
			t0_check[i][j] = ringQP.NewPoly()
			for jj := 0; jj < sizeB; jj++ {
				ringQP.NTT(B[j][jj], B[j][jj])
				ringQP.NTT(z[i][jj], z[i][jj])
				ringQP.MulCoeffsAndAdd(B[j][jj], z[i][jj], t0_check[i][j]) // Coefficient-wise multiplication
				ringQP.InvNTT(B[j][jj], B[j][jj])
				ringQP.InvNTT(z[i][jj], z[i][jj])
			}
			ringQP.InvNTT(t0_check[i][j], t0_check[i][j])

			ringQP.Sub(t0_check[i][j], w[i][j], t0_check[i][j])
			ringQP.Neg(t0_check[i][j], t0_check[i][j])

			ringQP.NTT(t[j], t[j])
			ringQP.NTT(t0_check[i][j], t0_check[i][j])
			ringQP.MulCoeffsAndAdd(tmpC, t[j], t0_check[i][j])
			ringQP.InvNTT(t[j], t[j])
			ringQP.InvNTT(t0_check[i][j], t0_check[i][j])

			// Check the B.z - w - sigma(c)t0 =?= 0
			temp.Zero()
			if t0_check[i][j].Equals(temp) == false {
				fmt.Printf("B.z != w + sigma(c)t0  for k=%v, j=%v \n", i, j)
				panic("WRONG CHECK -- B.z != w + sigma(c)t0")
			}
		}
	}
	fmt.Printf("[CHECK] B.z =?= w + sigma(c)t0  -- OK \n")

	// Check h_(0)=...=h_(k-1)=0
	for i := 0; i < k; i++ {
		if h.Coeffs[0][i] != 0 {
			panic("wrong h_0, ..., h_(k-1) != 0")
		} else {
			continue
		}
	}
	fmt.Printf("[CHECK] h_(0)=...=h_(k-1)=0     -- OK\n")

	// Create tau
	// Create iNTT(N.Atgamma ° t1) -- In Poly form
	v_prodk := make([][]*ring.Poly, k)
	for i := 0; i < k; i++ {
		v_prodk[i] = make([]*ring.Poly, m1)
		for j := 0; j < m1; j++ {
			v_prodk[i][j] = ringQP.NewPoly()
			ringQP.NTT(t[n+j], t[n+j])
			ringQP.MulCoeffs(Atgammak[i][j], t[n+j], v_prodk[i][j])
			ringQP.MulScalar(v_prodk[i][j], N, v_prodk[i][j])
			ringQP.InvNTT(t[n+j], t[n+j])
			ringQP.InvNTT(v_prodk[i][j], v_prodk[i][j])
		}
	}

	// Sum_v=0^m1 [iNTT(N.Atgamma ° t1)] - <u, gamma>
	for i := 0; i < k; i++ {
		for j := 1; j < m1; j++ {
			ringQP.Add(v_prodk[i][0], v_prodk[i][j], v_prodk[i][0])
		}
		ringQP.Add(v_prodk[i][0], cste[i], v_prodk[i][0])
	}

	// Create the target value T = Sum_v=0^k-1 \sig_v[ prod_k ]
	v_sumSig := make([]*ring.Poly, k)
	for i := 0; i < k; i++ {
		v_sumSig[i] = ringQP.NewPoly()

		// TODO sig is not a rotation!!!
		ringQP.Add(v_sumSig[i], v_prodk[i][0], v_sumSig[i])

	}

	// Create the value tau = Sum_v X^v T_v
	tau := ringQP.NewPoly()
	tau = v_sumSig[0].CopyNew()
	for i := 0; i < k; i++ {
		// TODO k>1
	}

	// Check tau equation
	// Create the scalar product v_scalarProdVm = <Atgamma_b1m[0], z1> TODO: m1>1; k>1
	v_scalarProdVm := make([][][][]*ring.Poly, k)
	for i := 0; i < k; i++ {
		v_scalarProdVm[i] = make([][][]*ring.Poly, k)
		for mu := 0; mu < k; mu++ {
			v_scalarProdVm[i][mu] = make([][]*ring.Poly, k)
			for nu := 0; nu < k; nu++ {
				v_scalarProdVm[i][mu][nu] = make([]*ring.Poly, m1)
				for j := 0; j < m1; j++ {
					v_scalarProdVm[i][mu][nu][j] = ringQP.NewPoly()
					for jj := 0; jj < sizeB; jj++ {
						ringQP.NTT(Atgamma_b1m[mu][j][jj], Atgamma_b1m[mu][j][jj])
						ringQP.NTT(z[(k+i-nu)%k][jj], z[(k+i-nu)%k][jj])
						ringQP.MulCoeffsAndAdd(Atgamma_b1m[mu][j][jj], z[(k+i-nu)%k][jj], v_scalarProdVm[i][mu][nu][j])
						ringQP.InvNTT(Atgamma_b1m[mu][j][jj], Atgamma_b1m[mu][j][jj])
						ringQP.InvNTT(z[(k+i-nu)%k][jj], z[(k+i-nu)%k][jj])
					}
					ringQP.InvNTT(v_scalarProdVm[i][mu][nu][j], v_scalarProdVm[i][mu][nu][j])
				}
			}
		}
	}

	// Sum the sum_mu sum_nu sum_j v_scalarProdVm[i][mu][nu][j]
	v_scalarProdV := make([]*ring.Poly, k)
	for i := 0; i < k; i++ {
		v_scalarProdV[i] = ringQP.NewPoly()

		for mu := 0; mu < k; mu++ {
			for nu := 0; nu < k; nu++ {
				for j := 0; j < m1; j++ {

					// TODO 1/k.X^mu

					// TODO sig^nu

					ringQP.Add(v_scalarProdV[i], v_scalarProdVm[i][mu][nu][j], v_scalarProdV[i])
				}
			}
		}
	}

	// Create the scalar product <b2, zi>
	scalarProd_b2z := make([]*ring.Poly, k)
	for i := 0; i < k; i++ {
		scalarProd_b2z[i] = ringQP.NewPoly()
		for j := 0; j < sizeB; j++ {
			ringQP.NTT(bVec[m1][j], bVec[m1][j])
			ringQP.NTT(z[i][j], z[i][j])
			ringQP.MulCoeffsAndAdd(bVec[m1][j], z[i][j], scalarProd_b2z[i])
			ringQP.InvNTT(bVec[m1][j], bVec[m1][j])
			ringQP.InvNTT(z[i][j], z[i][j])
		}
		ringQP.InvNTT(scalarProd_b2z[i], scalarProd_b2z[i])
	}

	// Create v_vi = v_scalarProdV[i] + <b2, zi>
	v_vVector := make([]*ring.Poly, k)
	for i := 0; i < k; i++ {
		v_vVector[i] = ringQP.NewPoly()
		ringQP.Add(v_scalarProdV[i], scalarProd_b2z[i], v_vVector[i])
	}

	// Create v_tauEquation = vi + sigma^i(c))(tau + t2 - h)
	v_tauEquation := make([]*ring.Poly, k)
	for i := 0; i < k; i++ {
		v_tauEquation[i] = ringQP.NewPoly()
		ringQP.Add(tau, t[n+m1], v_tauEquation[i])
		ringQP.Sub(v_tauEquation[i], h, v_tauEquation[i])

		// // TODO sig^i
		temp = c.CopyNew()

		ringQP.NTT(temp, temp)
		ringQP.NTT(v_tauEquation[i], v_tauEquation[i])
		ringQP.MulCoeffs(temp, v_tauEquation[i], v_tauEquation[i])
		ringQP.InvNTT(v_tauEquation[i], v_tauEquation[i])
		ringQP.InvNTT(temp, temp)
		ringQP.Add(v_tauEquation[i], vVector[i], v_tauEquation[i])
	}

	// Check v_tauEquation == v_vi
	for i := 0; i < k; i++ {
		if v_tauEquation[i].Equals(v_vVector[i]) == false {
			fmt.Printf("\n    [FAULT] v_tauEquation != v_vi  for k=%v \n", i)
			panic("[WRONG CHECK] -- v_tauEquation != v_vi")
		}
	}
	fmt.Printf("[CHECK] v_tauEquation =?= v_vi  -- OK \n")

}

func transpose(Ain [][]*ring.Poly) [][]*ring.Poly {
	lenX := len(Ain[0])
	lenY := len(Ain)
	Aout := make([][]*ring.Poly, lenX)
	for i := 0; i < lenX; i++ {
		Aout[i] = make([]*ring.Poly, lenY)
	}
	for i := 0; i < lenX; i++ {
		for j := 0; j < lenY; j++ {
			Aout[i][j] = Ain[j][i]
		}
	}
	return Aout
}

func main() {
	mainfunc()
}
