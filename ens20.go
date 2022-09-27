package main

import (
	"fmt"
	"github.com/ldsec/lattigo/v2/bfv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// Prover represents the prover party in the protocol.
//type Prover struct{}

// Verifier represents the verifier party in the protocol.
//type Verifier struct{}

// Protocol represents an instance of the protocol described in ENS20.
type Protocol struct {
	//Prover   Prover
	//Verifier Verifier

	uniformSampler  *ring.UniformSampler
	ternarySampler  *ring.TernarySampler
	gaussianSampler *ring.GaussianSampler
}

type PublicParams struct {
	N         int // dimension of Rq
	n         int // dimension for MSIS
	m1        int // dimension of the input in Rq
	m2        int // dimension of the input in Rq
	lbd       int // dimension for MLWE
	sizeB     int // dimension of bi
	size_t    int // dimension of the commitment
	size_bVec int // dimension of bVec vector of {bi}
	k         int // repeat rate in ENS20
	d1        int // delta1 bound on the mask // Change for appropriate value
}

func Execute(msg []*ring.Poly, publicParams PublicParams) (*Protocol, error) {
	// Initialize the ring parameters.
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParams, err := bfv.NewParametersFromLiteral(ringParamDef)
	if err != nil {
		return nil, fmt.Errorf("could not initialize the ring parameters: %s", err)
	}
	// Initialize the ring.
	ringQP := ringParams.RingQP()
	// Create the samplers.
	prng, err := utils.NewPRNG()
	if err != nil {
		return nil, fmt.Errorf("could not initialize the prng: %s", err)
	}
	uniformSampler := ring.NewUniformSampler(prng, ringQP)
	ternarySampler := ring.NewTernarySampler(prng, ringQP, 1.0/3.0, false)
	gaussianSampler := ring.NewGaussianSampler(prng, ringQP, ringParams.Sigma(), publicParams.d1)
	// Split the message over the polynomial space.
	mHat := make([]*ring.Poly, publicParams.m1)
	for i := 0; i < publicParams.m1; i++ {
		mHat[i] = ringQP.NewPoly()
		ringQP.InvNTT(msg[i], mHat[i])
	}

	// Construct the public matrix A
	// TODO: provide from the outside
	A := RandomMatrix(publicParams.m2, publicParams.m1, ringQP, uniformSampler)
	// Calculate Am = u
	u := MatrixVectorMul(A, msg, PolyMult{}, ringQP)
	// Sample a polynomial g s.t. g_0=...=g_(k-1)=0
	g := ringQP.NewPoly()
	uniformSampler.Read(g)
	for i := 0; i < publicParams.k; i++ {
		g.Coeffs[0][i] = 0
	}
	// Construct the matrix B.
	B := RandomMatrix(publicParams.n, publicParams.sizeB, ringQP, uniformSampler)
	// Construct the b_i vectors.
	bVecs := RandomMatrix(publicParams.size_bVec, publicParams.sizeB, ringQP, uniformSampler)
	// Construct the randomness r.
	r := RandomVector(publicParams.sizeB, ringQP, ternarySampler)
	// t_0 = Br.
	t_0 := MatrixVectorMul(B, r, PolyMult{}, ringQP)
	// t_1[i] = <b1[k], r> + mHat[k]  where mHat is the invNTT transform of msg
	t_1 := make([]*ring.Poly, len(mHat))
	for i := 0; i < len(t_1); i++ {
		t_1[i] = DotProduct(bVecs[i], r, PolyMult{}, ringQP)
		ringQP.Add(t_1[i], mHat[i], t_1[i])
	}
	// t2 = <b2, r> + g
	t_2 := DotProduct(bVecs[publicParams.m1], r, PolyMult{}, ringQP)
	ringQP.Add(t_2, g, t_2)
	// t = t0 || t1 || t2
	t := append(t_0, append(t_1, t_2)...)
	// Create the masks y, w
	y := RandomMatrix(publicParams.k, publicParams.sizeB, ringQP, gaussianSampler)
	w := make([][]*ring.Poly, publicParams.k)
	for i := 0; i < publicParams.k; i++ {
		w[i] = MatrixVectorMul(B, y[i], PolyMult{}, ringQP)
	}

	// Create the gamma challenges.
	gamma := RandomMatrix(publicParams.k, publicParams.m2, ringQP, uniformSampler)

	// Create h, the masked g
	At := MatrixTranspose(A)
	// Create the target value (At*gamma_k) in NTT form
	Atgammak_NTT := make([][]*ring.Poly, publicParams.k)
	for i := 0; i < publicParams.k; i++ {
		Atgammak_NTT[i] = make([]*ring.Poly, publicParams.m1)
		MatrixVectorMul(At, gamma[i], NTTMult{}, ringQP)
	}
	// Compute the inner product cste[i]=-<u, gamma[i]>  -- constant Poly
	temp := ringQP.NewPoly()
	cste := make([]*ring.Poly, publicParams.k)
	for i := 0; i < len(cste); i++ {
		cste[i] = DotProduct(u, gamma[i], NTTMult{}, ringQP)
		// --- ???: START
		for j := 0; j < ringParams.LogN(); j++ {
			ringQP.Shift(cste[i], 1<<j, temp)
			ringQP.Add(cste[i], temp, cste[i])
		}
		tmpInner := cste[i].Coeffs[0][0]
		cste[i].Zero()
		cste[i].Coeffs[0][0] = tmpInner
		// Negate the inner product
		ringQP.Neg(cste[i], cste[i])
		ringQP.Reduce(cste[i], cste[i])
		// --- ???: END
	}
	// Create  iNTT(N.Atgamma ° m[j]) = iNTT(N.Atgamma)* m^[j] -- in Poly form
	prodk := make([][]*ring.Poly, publicParams.k)
	for i := 0; i < publicParams.k; i++ {
		prodk[i] = MatrixVectorMul(Atgammak_NTT, msg, NTTMult{}, ringQP)
		for j := 0; j < publicParams.m1; j++ {
			ringQP.MulScalar(prodk[i][j], uint64(publicParams.N), prodk[i][j])
			ringQP.InvNTT(prodk[i][j], prodk[i][j])
		}
	}
	// Sum_v=0^m1 [iNTT(N.Atgamma ° m_v)] - <u, gamma> -- In Poly form
	for i := 0; i < publicParams.k; i++ {
		prodk[i][0] = VectorSum(prodk[i], ringQP)
		ringQP.Add(prodk[i][0], cste[i], prodk[i][0])
	}
	// Create the target value T = Sum_v=0^k-1 \sig_v[ prod_k ] -- In Poly form
	sumSig := make([]*ring.Poly, publicParams.k)
	for i := 0; i < publicParams.k; i++ {
		sumSig[i] = ringQP.NewPoly()
		// TODO sig is not a rotation!!!
		ringQP.Add(sumSig[i], prodk[i][0], sumSig[i])
	}
	// Create the value f = Sum_k X^k T[k]  -- In Poly form TODO
	fValue := ringQP.NewPoly()
	fValue = sumSig[0].CopyNew()
	for i := 0; i < publicParams.k; i++ {
		// TODO
	}
	h := ringQP.NewPoly()
	ringQP.Add(g, fValue, h)

	// Create v_i
	// Create the value Atgamma_b1m[i][j][k] = iNTT(N.Atgamma[i][j] ° b1[k]) for k in sizeB  -- in Poly form
	Atgamma_b1m := make([][][]*ring.Poly, publicParams.k)
	bVecs_NTT := ApplyMatrix(bVecs, ToNTT{}, ringQP)
	for i := 0; i < len(Atgamma_b1m); i++ {
		Atgamma_b1m[i] = make([][]*ring.Poly, publicParams.m1)
		for j := 0; j < len(Atgamma_b1m[0]); j++ {
			Atgamma_b1m[i][j] = make([]*ring.Poly, publicParams.sizeB)
			for k := 0; k < len(Atgamma_b1m[0][0]); k++ {
				Atgamma_b1m[i][j][k] = NTTMult{}.Apply(Atgammak_NTT[i][j], bVecs_NTT[j][k], ringQP)
				ringQP.MulScalar(Atgamma_b1m[i][j][k], uint64(publicParams.N), Atgamma_b1m[i][j][k])
				ringQP.InvNTT(Atgamma_b1m[i][j][k], Atgamma_b1m[i][j][k])
			}
		}
	}
	// Create the scalar product scalarProdVm = <Atgamma_b1m[0]
	scalarProdVm := NewMultiArray([]int{publicParams.k, publicParams.k, publicParams.k, publicParams.m1, publicParams.sizeB}, ringQP)
	scalarProdVm.Map(func(el *ring.Poly, coords []int) {
		i, mu, nu, j, k := coords[0], coords[1], coords[2], coords[3], coords[4]
		yP := (publicParams.k + i - nu) % publicParams.k
		DotProduct()
		y[yP]
		ringQP.NTT(Atgamma_b1m[mu][j][k], Atgamma_b1m[mu][j][k])
		ringQP.NTT(y[yP][k], y[yP][k])
		ringQP.MulCoeffsAndAdd(Atgamma_b1m[mu][j][k], y[yP][k], scalarProdVm[i][mu][nu][j])
		ringQP.InvNTT(Atgamma_b1m[mu][j][k], Atgamma_b1m[mu][j][k])
		ringQP.InvNTT(y[yP][k], y[yP][k])
	})
	scalarProdVm := make([][][][]*ring.Poly, publicParams.k)
	for i := 0; i < publicParams.k; i++ {
		scalarProdVm[i] = make([][][]*ring.Poly, publicParams.k)
		for mu := 0; mu < publicParams.k; mu++ {
			scalarProdVm[i][mu] = make([][]*ring.Poly, publicParams.k)
			for nu := 0; nu < publicParams.k; nu++ {
				scalarProdVm[i][mu][nu] = make([]*ring.Poly, publicParams.m1)
				for j := 0; j < publicParams.m1; j++ {
					scalarProdVm[i][mu][nu][j] = ringQP.NewPoly()
					for k := 0; k < publicParams.sizeB; k++ {
						yP := (publicParams.k + i - nu) % publicParams.k
						ringQP.NTT(Atgamma_b1m[mu][j][k], Atgamma_b1m[mu][j][k])
						ringQP.NTT(y[yP][k], y[yP][k])
						ringQP.MulCoeffsAndAdd(Atgamma_b1m[mu][j][k], y[yP][k], scalarProdVm[i][mu][nu][j])
						ringQP.InvNTT(Atgamma_b1m[mu][j][k], Atgamma_b1m[mu][j][k])
						ringQP.InvNTT(y[yP][k], y[yP][k])
					}
					ringQP.InvNTT(scalarProdVm[i][mu][nu][j], scalarProdVm[i][mu][nu][j])
				}
			}
		}
	}

	// Sum the sum_mu sum_nu sum_j scalarProdVm[i][mu][nu][j]
	scalarProdV := make([]*ring.Poly, publicParams.k)
	for i := 0; i < publicParams.k; i++ {
		scalarProdV[i] = ringQP.NewPoly()
		for mu := 0; mu < publicParams.k; mu++ {
			for nu := 0; nu < publicParams.k; nu++ {
				for j := 0; j < publicParams.m1; j++ {
					// TODO 1/k.X^mu
					// TODO sig^nu
					ringQP.Add(scalarProdV[i], scalarProdVm[i][mu][nu][j], scalarProdV[i])
				}
			}
		}
	}

	// Create the scalar product <b2, yi>
	scalarProd_b2y := make([]*ring.Poly, publicParams.k)
	for i := 0; i < len(scalarProd_b2y); i++ {
		scalarProd_b2y[i] = DotProduct(bVecs[publicParams.m1], y[i], PolyMult{}, ringQP)
	}

	// Create vi = scalarProdV + <b2, yi>
	vVector := VectorAdd(scalarProdV, scalarProd_b2y, ringQP)

	// Get challenge c in C: {-1, 0, 1}^N
	c := ringQP.NewPoly()
	ternarySampler.Read(c)

	// Create the z_i openings
	var z = make([][]*ring.Poly, publicParams.k)
	for i := 0; i < publicParams.k; i++ {
		z[i] = make([]*ring.Poly, publicParams.sizeB)

		for j := 0; j < publicParams.sizeB; j++ {
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

	return nil, nil
}
