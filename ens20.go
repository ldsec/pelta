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

func Execute(msg MultiArray, publicParams PublicParams) (*Protocol, error) {
	// Initialize the ring parameters.
	ringParamDef := bfv.PN13QP218
	ringParamDef.T = 0x3ee0001
	ringParams, err := bfv.NewParametersFromLiteral(ringParamDef)
	if err != nil {
		return nil, fmt.Errorf("could not initialize the ring parameters: %s", err)
	}
	// Initialize the ring.
	baseRing := ringParams.RingQP()
	// Create the samplers.
	prng, err := utils.NewPRNG()
	if err != nil {
		return nil, fmt.Errorf("could not initialize the prng: %s", err)
	}
	uniformSampler := ring.NewUniformSampler(prng, baseRing)
	ternarySampler := ring.NewTernarySampler(prng, baseRing, 1.0/3.0, false)
	gaussianSampler := ring.NewGaussianSampler(prng, baseRing, ringParams.Sigma(), publicParams.d1)
	// Split the message over the polynomial space.
	mHat := NewVectorFromDimensions(publicParams.m1, baseRing)
	mHat.ForEach(func(el *ring.Poly, coords []int) {
		baseRing.InvNTT(el, el)
	})
	// Construct the public matrix A
	// TODO: provide from the outside
	A := NewRandomMatrix(publicParams.m2, publicParams.m1, baseRing, uniformSampler)
	// Calculate Am = u
	u := MatrixVectorMul(A, msg, PolyMult{}, baseRing)
	// Sample a polynomial g s.t. g_0=...=g_(k-1)=0
	g := baseRing.NewPoly()
	uniformSampler.Read(g)
	for i := 0; i < publicParams.k; i++ {
		g.Coeffs[0][i] = 0
	}
	// Construct the matrix B.
	B := NewRandomMatrix(publicParams.n, publicParams.sizeB, baseRing, uniformSampler)
	// Construct the b_i vectors.
	bVecs := NewRandomMatrix(publicParams.size_bVec, publicParams.sizeB, baseRing, uniformSampler)
	// Construct the randomness r.
	r := NewRandomVector(publicParams.sizeB, baseRing, ternarySampler)
	// t_0 = Br.
	t_0 := MatrixVectorMul(B, r, PolyMult{}, baseRing)
	// t_1[i] = <b1[k], r> + mHat[k]  where mHat is the invNTT transform of msg
	t_1 := NewVectorFromDimensions(mHat.Length(), baseRing)
	t_1.ForEach(func(el *ring.Poly, coords []int) {
		i := coords[0]
		tmp := DotProduct(bVecs.MatrixRowSlice(i), r, PolyMult{}, baseRing)
		baseRing.Add(tmp, mHat.ElementAtIndex(i), el)
	})
	// t2 = <b2, r> + g
	t_2 := DotProduct(bVecs.MatrixRowSlice(publicParams.m1), r, PolyMult{}, baseRing)
	baseRing.Add(t_2, g, t_2)
	// t = t0 || t1 || t2
	t := NewVectorFromSlice(append(t_0.Array, append(t_1.Array, t_2)...))
	// Create the masks y, w
	y := NewRandomMatrix(publicParams.k, publicParams.sizeB, baseRing, gaussianSampler)
	w := NewMatrixFromDimensions(publicParams.n, publicParams.k, baseRing)
	for i := 0; i < publicParams.k; i++ {
		w.MatrixSetRow(i, MatrixVectorMul(B, y.MatrixRowSlice(i), PolyMult{}, baseRing))
	}
	// Create the gamma challenges.
	gamma := NewRandomMatrix(publicParams.k, publicParams.m2, baseRing, uniformSampler)
	// Create h, the masked g
	At := A.MatrixTranspose()
	// Create the target value (At*gamma_k) in NTT form
	Atgammak_NTT := NewMatrixFromDimensions(publicParams.k, publicParams.m1, baseRing)
	for i := 0; i < publicParams.k; i++ {
		Atgammak_NTT.MatrixSetRow(i, MatrixVectorMul(At, gamma.MatrixRowSlice(i), NTTMult{}, baseRing))
	}
	// Compute the inner product cste[i]=-<u, gamma[i]>  -- constant Poly
	temp := baseRing.NewPoly()
	cste := NewVectorFromDimensions(publicParams.k, baseRing)
	cste.ForEach(func(el *ring.Poly, coords []int) {
		i := coords[0]
		dt := DotProduct(u, gamma.MatrixRowSlice(i), NTTMult{}, baseRing)
		// --- ???: START
		for j := 0; j < ringParams.LogN(); j++ {
			baseRing.Shift(dt, 1<<j, temp)
			baseRing.Add(dt, temp, el)
		}
		tmpInner := dt.Coeffs[0][0]
		dt.Zero()
		dt.Coeffs[0][0] = tmpInner
		// Negate the inner product
		baseRing.Neg(dt, dt)
		baseRing.Reduce(dt, dt)
		// --- ???: END
		el.Copy(dt)
	})
	// Create  iNTT(N.Atgamma ° m[j]) = iNTT(N.Atgamma)* m^[j] -- in Poly form
	prodk := NewMatrixFromDimensions(publicParams.k, publicParams.k, baseRing)
	for i := 0; i < publicParams.k; i++ {
		prodk.MatrixSetRow(i, MatrixVectorMul(Atgammak_NTT, msg, NTTMult{}, baseRing))
	}
	prodk.ForEach(func(el *ring.Poly, _ []int) {
		baseRing.MulScalar(el, uint64(publicParams.N), el)
		baseRing.InvNTT(el, el)
	})
	// Sum_v=0^m1 [iNTT(N.Atgamma ° m_v)] - <u, gamma> -- In Poly form
	for i := 0; i < publicParams.k; i++ {
		prodk[i][0] = VectorSum(prodk[i], baseRing)
		baseRing.Add(prodk[i][0], cste[i], prodk[i][0])
	}
	// Create the target value T = Sum_v=0^k-1 \sig_v[ prod_k ] -- In Poly form
	// TODO: use the new sigma function
	sumSig := NewVectorFromDimensions(publicParams.k, baseRing)
	// Create the value f = Sum_k X^k T[k]  -- In Poly form TODO
	fValue := baseRing.NewPoly()
	fValue = sumSig.Array[0].CopyNew()
	for i := 0; i < publicParams.k; i++ {
		// TODO
	}
	h := baseRing.NewPoly()
	baseRing.Add(g, fValue, h)

	// Create v_i
	// Create the value Atgamma_b1m[i][j][k] = iNTT(N.Atgamma[i][j] ° b1[k]) for k in sizeB  -- in Poly form
	Atgamma_b1m := NewMultiArray([]int{publicParams.k, publicParams.m1, publicParams.sizeB}, baseRing)
	// Convert bVecs into NTT form.
	bVecs.ForEach(func(el *ring.Poly, _ []int) { baseRing.NTT(el, el) })
	Atgamma_b1m.ForEach(func(el *ring.Poly, coords []int) {
		i, j, k := coords[0], coords[1], coords[2]
		NTTMult{}.Apply(Atgammak_NTT.MatrixElement(i, j), bVecs.MatrixElement(j, k), el, baseRing)
		baseRing.MulScalar(el, uint64(publicParams.N), el)
		baseRing.InvNTT(el, el)
	})
	// Convert bVecs back into poly form.
	bVecs.ForEach(func(el *ring.Poly, _ []int) { baseRing.InvNTT(el, el) })
	// Create the scalar product scalarProdVm = <Atgamma_b1m[0]
	scalarProdVm := NewMultiArray([]int{publicParams.k, publicParams.k, publicParams.k, publicParams.m1}, baseRing)
	scalarProdVm.ForEach(func(el *ring.Poly, coords []int) {
		i, mu, nu, j := coords[0], coords[1], coords[2], coords[3]
		for k := 0; k < publicParams.sizeB; k++ {
			// TODO compress into NTT dot product
			yP := y.MatrixRowSlice((publicParams.k + i - nu) % publicParams.k)
			target := Atgamma_b1m.ElementAtCoords([]int{mu, j, k})
			baseRing.NTT(target, target)
			baseRing.NTT(yP.ElementAtIndex(k), yP.ElementAtIndex(k))
			baseRing.MulCoeffsAndAdd(target, yP.ElementAtIndex(k), el)
			baseRing.InvNTT(target, target)
			baseRing.InvNTT(yP.ElementAtIndex(k), yP.ElementAtIndex(k))
		}
	})
	scalarProdV := NewVectorFromDimensions(publicParams.k, baseRing)
	scalarProdV.ForEach(func(el1 *ring.Poly, _ []int) {
		scalarProdVm.ForEach(func(el2 *ring.Poly, _ []int) {
			// TODO 1/k.X^mu
			// TODO sig^nu
			baseRing.Add(el1, el2, el1)
		})
	})

	// Create the scalar product <b2, yi>
	scalarProd_b2y := NewVectorFromDimensions(publicParams.k, baseRing)
	for i := 0; i < scalarProd_b2y.Length(); i++ {
		scalarProd_b2y.SetElementAtIndex(i,
			DotProduct(bVecs.MatrixRowSlice(publicParams.m1), y.MatrixRowSlice(i), PolyMult{}, baseRing))
	}

	// Create vi = scalarProdV + <b2, yi>
	vVector := VectorAdd(scalarProdV, scalarProd_b2y, baseRing)

	// Get challenge c in C: {-1, 0, 1}^N
	c := baseRing.NewPoly()
	ternarySampler.Read(c)

	// Create the z_i openings
	var z = make([][]*ring.Poly, publicParams.k)
	for i := 0; i < publicParams.k; i++ {
		z[i] = make([]*ring.Poly, publicParams.sizeB)

		for j := 0; j < publicParams.sizeB; j++ {
			z[i][j] = baseRing.NewPoly()
			factor := baseRing.NewPoly()
			//baseRing.Shift // TODO sigma not a shift
			factor = c.CopyNew()
			baseRing.NTT(r[j], r[j])
			baseRing.NTT(factor, factor)
			baseRing.MulCoeffs(factor, r[j], z[i][j])
			baseRing.InvNTT(factor, factor)
			baseRing.InvNTT(r[j], r[j])
			baseRing.InvNTT(z[i][j], z[i][j])

			baseRing.Add(z[i][j], y[i][j], z[i][j])
		}
	}

	return nil, nil
}
