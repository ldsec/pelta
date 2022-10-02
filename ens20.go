package main

import (
	"fmt"

	"github.com/ldsec/codeBase/commitment/math"
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
	n         int // dimension for M-SIS
	m1        int // dimension of the input in Rq
	m2        int // dimension of the input in Rq
	lbd       int // dimension for M-LWE
	sizeB     int // dimension of bi
	size_t    int // dimension of the commitment
	size_bVec int // dimension of bVec vector of {bi}
	k         int // repeat rate in ENS20
	d1        int // delta1 bound on the mask // Change for appropriate value
}

type PublicParams2 struct {
	d      int // deg(X^d + 1), a power of two
	q      int // Rational prime mod
	m      int // # rows
	n      int // # cols
	k      int // Repetition rate
	delta1 int // Width of the uniform distribution
	lambda int // M-LWE dimension
	kappa  int // M-SIS dimension
}

func Execute2(s math.Vector, params PublicParams2) (*Protocol, error) {

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
	gaussianSampler := ring.NewGaussianSampler(prng, baseRing, ringParams.Sigma(), params.delta1)

	// Inputs
	A := NewRandomIntegerMatrix(params.m, params.n, params.q)
	u := A.DeepCopy().AsMatrix().MulVec(s)
	B0 := NewRandomPolynomialMatrix(params.kappa, params.lambda+params.kappa+params.n/params.d+3, baseRing, uniformSampler)

	return nil, nil
}

func Execute(msg math.Vector, publicParams PublicParams) (*Protocol, error) {
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
	mHat := math.NewVectorFromSize(publicParams.m1).AsPolyArray().InvNTT().AsVector()
	// Construct the public matrix A
	// TODO: provide from the outside
	A := NewRandomPolynomialMatrix(publicParams.m2, publicParams.m1, baseRing, uniformSampler)
	// Calculate Am = u
	u := A.MulVec(msg)
	// Sample a polynomial g s.t. g_0=...=g_(k-1)=0
	g := math.NewPolynomial(baseRing)
	uniformSampler.Read(g.Ref)
	for i := 0; i < publicParams.k; i++ {
		g.Ref.Coeffs[0][i] = 0
	}
	// Construct the matrix B.
	B := NewRandomPolynomialMatrix(publicParams.n, publicParams.sizeB, baseRing, uniformSampler)
	// Construct the b_i vectors.
	bVecs := NewRandomPolynomialMatrix(publicParams.size_bVec, publicParams.sizeB, baseRing, uniformSampler)
	// Construct the randomness r.
	r := NewRandomPolynomialVector(publicParams.sizeB, baseRing, ternarySampler)
	// t_0 = Br.
	t_0 := B.MulVec(r)
	// t_1[i] = <b1[k], r> + mHat[k]  where mHat is the invNTT transform of msg
	t_1 := math.NewVectorFromSize(mHat.Length()).Populate(
		func(i int) math.RingElement {
			return bVecs.Row(i).DeepCopy().AsPolyArray().
				NTT().AsVector().
				DotProduct(r).(math.Polynomial).
				InvNTT().
				Add(mHat.ElementAtIndex(i))
		})
	// t2 = <b2, r> + g
	t_2 := bVecs.Row(publicParams.m1).DeepCopy().AsVector().DotProduct(r).Add(g)
	// t = t0 || t1 || t2
	t := math.NewVectorFromSlice(append(t_0.Array, append(t_1.Array, t_2)...))
	// Create the masks y, w
	y := NewRandomPolynomialMatrix(publicParams.k, publicParams.sizeB, baseRing, gaussianSampler)
	w := math.NewMatrixFromDimensions(publicParams.k, publicParams.n).PopulateRows(func(i int) math.Vector {
		return B.MulVec(y.Row(i))
	})
	// Create the gamma challenges.
	gamma := NewRandomPolynomialMatrix(publicParams.k, publicParams.m2, baseRing, uniformSampler)
	// Create h, the masked g
	At := A.DeepCopy().AsMatrix().Transpose()
	// Create the target value (At*gamma_k) in NTT form
	Atgammak_NTT := math.NewMatrixFromDimensions(publicParams.k, publicParams.m1).PopulateRows(func(i int) math.Vector {
		return At.DeepCopy().AsMatrix().MulVec(gamma.Row(i))
	})
	// Compute the inner product cste[i]=-<u, gamma[i]>  -- constant Poly
	temp := baseRing.NewPoly()
	cste := math.NewVectorFromSize(publicParams.k)
	cste.ForEach(func(el math.RingElement, i int) {
		dt := u.DotProduct(gamma.Row(i)).(math.Polynomial)
		// --- ???: START
		for j := 0; j < ringParams.LogN(); j++ {
			baseRing.Shift(dt.Ref, 1<<j, temp)
			baseRing.Add(dt.Ref, temp, el.(math.Polynomial).Ref)
		}
		tmpInner := dt.Ref.Coeffs[0][0]
		dt.Ref.Zero()
		dt.Ref.Coeffs[0][0] = tmpInner
		// Negate the inner product
		baseRing.Neg(dt.Ref, dt.Ref)
		baseRing.Reduce(dt.Ref, dt.Ref)
		// --- ???: END
		el.(math.Polynomial).Ref.Copy(dt.Ref)
	})
	// Create  iNTT(N.Atgamma ° m[j]) = iNTT(N.Atgamma)* m^[j] -- in Poly form
	prodk := math.NewMatrixFromDimensions(publicParams.k, publicParams.k).PopulateRows(func(i int) math.Vector {
		return Atgammak_NTT.DeepCopy().AsMatrix().MulVec(msg).AsPolyArray().Scale(uint64(publicParams.N)).InvNTT().AsVector()
	})
	// Sum_v=0^m1 [iNTT(N.Atgamma ° m_v)] - <u, gamma> -- In Poly form
	for i := 0; i < publicParams.k; i++ {
		prodk.SetElement(i, 0, prodk.Row(i).Sum().Add(cste.ElementAtIndex(i)))
	}
	// Create the target value T = Sum_v=0^k-1 \sig_v[ prod_k ] -- In Poly form
	// TODO: use the new sigma function
	sumSig := math.NewVectorFromSize(publicParams.k)
	// Create the value f = Sum_k X^k T[k]  -- In Poly form TODO
	fValue := sumSig.Array[0].Copy()
	for i := 0; i < publicParams.k; i++ {
		// TODO
	}
	h := g.Copy().Add(fValue)
	// Create v_i
	// Convert bVecs into NTT form.
	bVecs = bVecs.AsPolyArray().NTT().AsMatrix()
	// Create the value Atgamma_b1m[i][j][k] = iNTT(N.Atgamma[i][j] ° b1[k]) for k in sizeB  -- in Poly form
	Atgamma_b1m := math.NewMultiArray([]int{publicParams.k, publicParams.m1, publicParams.sizeB})
	Atgamma_b1m.Map(func(el math.Polynomial, coords []int) math.Polynomial {
		i, j, k := coords[0], coords[1], coords[2]
		return Atgammak_NTT.Element(i, j).Mul(bVecs.Element(j, k)).(math.Polynomial).Scale(uint64(publicParams.N)).InvNTT()
	})
	// Convert bVecs back into poly form.
	bVecs = bVecs.AsPolyArray().InvNTT().AsMatrix()
	// Create the scalar product scalarProdVm = <Atgamma_b1m[0]
	scalarProdVm := math.NewMultiArray([]int{publicParams.k, publicParams.k, publicParams.k, publicParams.m1})
	scalarProdVm.ForEach(func(el math.RingElement, coords []int) {
		i, mu, nu, j := coords[0], coords[1], coords[2], coords[3]
		for k := 0; k < publicParams.sizeB; k++ {
			// TODO compress into NTT dot product
			yP := y.Row((publicParams.k + i - nu) % publicParams.k)
			target := Atgamma_b1m.ElementAtCoords([]int{mu, j, k}).(math.Polynomial).NTT()
			yP.ElementAtIndex(k).(math.Polynomial).NTT()
			baseRing.MulCoeffsAndAdd(target.Ref, yP.ElementAtIndex(k).(math.Polynomial).Ref, el.(math.Polynomial).Ref)
			target.InvNTT()
			yP.ElementAtIndex(k).(math.Polynomial).InvNTT()
		}
	})
	scalarProdV := math.NewVectorFromSize(publicParams.k)
	scalarProdV.ForEach(func(el1 math.RingElement, _ int) {
		scalarProdVm.ForEach(func(el2 math.RingElement, _ []int) {
			// TODO 1/k.X^mu
			// TODO sig^nu
			baseRing.Add(el1.(math.Polynomial).Ref, el2.(math.Polynomial).Ref, el1.(math.Polynomial).Ref)
		})
	})

	// Create the scalar product <b2, yi>
	scalarProd_b2y := math.NewVectorFromSize(publicParams.k).Populate(func(i int) math.RingElement {
		return bVecs.Row(publicParams.m1).DeepCopy().AsVector().DotProduct(y.Row(i))
	})
	// Create vi = scalarProdV + <b2, yi>
	vVector := scalarProdV.Add(scalarProd_b2y.MultiArray)

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
			baseRing.NTT(r.Array[j].(math.Polynomial).Ref, r.Array[j].(math.Polynomial).Ref)
			baseRing.NTT(factor, factor)
			baseRing.MulCoeffs(factor, r.Array[j].(math.Polynomial).Ref, z[i][j])
			baseRing.InvNTT(factor, factor)
			baseRing.InvNTT(r.Array[j].(math.Polynomial).Ref, r.Array[j].(math.Polynomial).Ref)
			baseRing.InvNTT(z[i][j], z[i][j])

			baseRing.Add(z[i][j], y.Element(i, j).(math.Polynomial).Ref, z[i][j])
		}
	}

	return nil, nil
}
