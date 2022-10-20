package main

//
//import (
//	"fmt"
//	"github.com/ldsec/codeBase/commitment/ens20"
//	"github.com/ldsec/codeBase/commitment/math"
//	"github.com/tuneinsight/lattigo/v4/bfv"
//	"github.com/tuneinsight/lattigo/v4/ring"
//	"github.com/tuneinsight/lattigo/v4/utils"
//)
//
//// Prover represents the prover party in the protocol.
////type Prover struct{}
//
//// Verifier represents the verifier party in the protocol.
////type Verifier struct{}
//
//// Protocol represents an instance of the protocol described in ENS20.
//type Protocol struct {
//	//Prover   Prover
//	//Verifier Verifier
//
//	uniformSampler  *ring.UniformSampler
//	ternarySampler  *ring.TernarySampler
//	gaussianSampler *ring.GaussianSampler
//}
//
//// ExecuteCubic executes the protocol for proving the knowledge of a ternary vector s as a solution to As = u.
//func ExecuteCubic(s math.Vector, settings Settings) (bool, error) {
//	// Initialize the ring parameters.
//	ringParamDef := bfv.PN13QP218
//	ringParamDef.T = 0x3ee0001
//	ringParams, err := bfv.NewParametersFromLiteral(ringParamDef)
//	if err != nil {
//		return false, fmt.Errorf("could not initialize the ring parameters: %s", err)
//	}
//	// The automorphism generator.
//	sig := math.NewAutomorphism(int64(settings.d), int64(settings.k)) //uint64(2*N/k+1)
//	// Initialize the ring.
//	baseRing := ringParams.RingQP().RingQ
//	// Create the samplers.
//	prng, err := utils.NewPRNG()
//	if err != nil {
//		return false, fmt.Errorf("could not initialize the prng: %s", err)
//	}
//	uniformSampler := ring.NewUniformSampler(prng, baseRing)
//	ternarySampler := ring.NewTernarySampler(prng, baseRing, 1.0/3.0, false)
//	gaussianSampler := ring.NewGaussianSampler(prng, baseRing, ringParams.Sigma(), settings.delta1)
//	// Inputs
//	numSplits := settings.n / settings.d
//	// Split the message into polynomial space.
//	sHat := ens20.SplitInvNTT(s, numSplits, settings.d, baseRing)
//	A := ens20.NewRandomIntegerMatrix(settings.m, settings.n, settings.q)
//	// As = u
//	u := A.Copy().AsMatrix().MulVec(s)
//	bSize := numSplits + 3
//	B0 := ens20.NewRandomPolynomialMatrix(settings.kappa, settings.lambda+settings.kappa+bSize, baseRing, uniformSampler)
//	b := ens20.NewRandomPolynomialMatrix(bSize, B0.Cols(), baseRing, uniformSampler)
//
//	// Sample a polynomial g s.t. g_0=...=g_(k-1)=0
//	g := math.NewZeroPolynomial(baseRing)
//	uniformSampler.Read(g.Ref)
//	for i := 0; i < settings.k; i++ {
//		g.SetCoefficient(i, 0)
//	}
//	// Sample the randomness.
//	r := ens20.NewRandomPolynomialVector(B0.Cols()*settings.d, baseRing, ternarySampler)
//	// Compute the commitments.
//	t0 := B0.Copy().AsMatrix().MulVec(r)
//	t := math.NewVectorFromSize(bSize).Populate(func(i int) math.RingElement {
//		if i > numSplits {
//			return math.NewZeroPolynomial(baseRing)
//		}
//		res := b.Row(i).Dot(r)
//		if i == numSplits {
//			return res.Add(g)
//		} else {
//			return res.Add(sHat.Element(i))
//		}
//	})
//
//	// Create the masks.
//	y := ens20.NewRandomPolynomialMatrix(settings.k, r.Length(), baseRing, gaussianSampler)
//	w := math.NewMatrixFromDimensions(settings.k, r.Length()).PopulateRows(func(i int) math.Vector {
//		return B0.Copy().AsMatrix().MulVec(y.Row(i))
//	})
//
//	// Verifier multiplicant generation.
//	alpha := ens20.NewRandomPolynomialVector(settings.k*numSplits, baseRing, uniformSampler)
//	gamma := ens20.NewRandomIntegerMatrix(settings.k, settings.m, settings.q)
//
//	// Prover further set up.
//	sum1 := math.NewMatrixFromDimensions(settings.k, numSplits).
//		Populate(func(i int, j int) math.RingElement {
//			tmp := sig.Permute(-1,
//				sHat.Element(j).Copy().(math.Polynomial).
//					Scale(uint64(3)).
//					Mul(b.Row(j).Dot(y.Row(i)).Pow(2)).(math.Polynomial))
//			return alpha.Element(i*numSplits + j).Copy().Mul(tmp)
//		}).Sum().Neg()
//	sum2 := math.NewMatrixFromDimensions(settings.k, numSplits).
//		Populate(func(i int, j int) math.RingElement {
//			tmp := sig.Permute(-1,
//				sHat.Element(j).Copy().(math.Polynomial).
//					Pow(2).(math.Polynomial).
//					Scale(3).
//					Add(math.NewZeroPolynomial(baseRing).One().Neg()).
//					Mul(b.Row(j).Dot(y.Row(i)).Pow(2)).(math.Polynomial))
//			return alpha.Element(i*numSplits + j).Copy().Mul(tmp)
//		}).Sum()
//	sum3 := math.NewMatrixFromDimensions(settings.k, numSplits).
//		Populate(func(i int, j int) math.RingElement {
//			tmp := sig.Permute(-1,
//				b.Row(j).
//					Dot(y.Row(i)).
//					Pow(3).(math.Polynomial))
//			return alpha.Element(i*numSplits + j).Copy().Mul(tmp)
//		}).Sum()
//	t.SetElementAtIndex(numSplits+1, b.Row(numSplits+1).Copy().AsVec().Dot(r).Add(b.Row(numSplits+2).Copy().AsVec().Dot(y.Row(0))).Add(sum1))
//	t.SetElementAtIndex(numSplits+2, b.Row(numSplits+2).Copy().AsVec().Dot(r).Add(sum2))
//	v := b.Row(numSplits + 1).Copy().AsVec().Dot(y.Row(0)).Add(sum3)
//	At := A.Copy().AsMatrix().Transpose()
//	psi := math.NewMatrixFromDimensions(settings.k, numSplits).PopulateRows(
//		func(mu int) math.Vector {
//			tmp := At.Copy().AsMatrix().MulVec(gamma.Row(mu))
//			return ens20.SplitInvNTT(tmp, numSplits, settings.d, baseRing)
//		})
//
//	invk := math.NewModInt(int64(settings.k), settings.q).Inv()
//	gMask := math.NewVectorFromSize(settings.k).Populate(
//		func(mu int) math.RingElement {
//			tmp := math.NewVectorFromSize(settings.k).Populate(
//				func(v int) math.RingElement {
//					dec := math.NewZeroPolynomial(baseRing).One().(math.Polynomial).Scale(u.Copy().AsVec().Dot(gamma.Row(mu)).(*math.ModInt).Uint64())
//					return sig.Permute(int64(v),
//						math.NewVectorFromSize(numSplits).Populate(
//							func(j int) math.RingElement {
//								return psi.Element(mu, j).Copy().Mul(sHat.Element(j)).(math.Polynomial).Scale(uint64(settings.d))
//							}).
//							Sum().
//							Add(dec.Neg()).(math.Polynomial))
//				}).Sum().(math.Polynomial)
//			return tmp.RRot(mu).Mul(invk)
//		}).Sum()
//	h := g.Add(gMask)
//	vp := math.NewVectorFromSize(settings.k).Populate(
//		func(i int) math.RingElement {
//			return math.NewVectorFromSize(settings.k).Populate(
//				func(mu int) math.RingElement {
//					tmp2 := math.NewMatrixFromDimensions(settings.k, numSplits).Populate(
//						func(v, j int) math.RingElement {
//							return sig.Permute(int64(v),
//								b.Row(j).
//									Mul(psi.Element(mu, j)).
//									Scale(uint64(settings.d)).AsVec().
//									Dot(y.Row(i-v)).(math.Polynomial))
//						}).Sum().(math.Polynomial)
//					return tmp2.RRot(mu).Mul(invk)
//				}).
//				Sum().
//				Add(b.Row(numSplits).Dot(y.Row(i)))
//		})
//
//	// Verifier challenge generation.
//	c := ens20.NewRandomPolynomial(baseRing, uniformSampler)
//
//	// Masked openings.
//	z := math.NewMatrixFromDimensions(settings.k, r.Length()).PopulateRows(
//		func(i int) math.Vector {
//			sigc := sig.Permute(int64(i), c.Copy().(math.Polynomial))
//			return y.Row(i).Copy().Add(r.Mul(sigc)).AsVec()
//		})
//
//	return VerifyCubic(settings, numSplits, t0, t, alpha, vp, u, A, B0, b, w, gamma, z, h.(math.Polynomial), v.(math.Polynomial), c, sig, baseRing), nil
//}
//
//func VerifyCubic(params Settings, numSplits int, t0, t, alpha, vp, u math.Vector, A, B0, b, w, gamma, z math.Matrix, h, v, c math.Polynomial, sig math.Automorphism, baseRing *ring.Ring) bool {
//	maskedOpeningTest := z.AllRows(
//		func(zi math.Vector, i int) bool {
//			rhs := w.Row(i).Copy().Add(t0.Mul(sig.Permute(int64(i), c)))
//			return zi.AsCoeffs().L2Norm() < params.beta && B0.MulVec(zi).Eq(rhs)
//		})
//	if !maskedOpeningTest {
//		return false
//	}
//	// Constructing f
//	f := math.NewMatrixFromDimensions(params.k, numSplits).Populate(
//		func(i int, j int) math.RingElement {
//			return b.Row(j).Dot(z.Row(i)).Add(t.Element(j).Copy().Mul(sig.Permute(int64(i), c)).Neg())
//		})
//	f2 := b.Row(numSplits + 2).Dot(z.Row(0)).Add(c.Mul(t.Element(numSplits + 2)).Neg())
//	f3 := b.Row(numSplits + 3).Dot(z.Row(0)).Add(c.Mul(t.Element(numSplits + 3)).Neg())
//	vTest := math.NewMatrixFromDimensions(params.k, numSplits).Populate(
//		func(i int, j int) math.RingElement {
//			p1 := f.Element(i, j).Copy()
//			p2 := f.Element(i, j).Copy().Add(sig.Permute(int64(i), c))
//			p3 := f.Element(i, j).Copy().Add(sig.Permute(int64(i), c).Neg())
//			return alpha.Element(numSplits*i + j).Mul(sig.Permute(int64(-i), p1.Mul(p2).Mul(p3).(math.Polynomial)))
//		}).Sum().
//		Add(f2).
//		Add(c.Mul(f3)).
//		Eq(v)
//	if !vTest {
//		return false
//	}
//	hTest := math.NewVectorFromSize(params.k).All(
//		func(_ math.RingElement, i int) bool {
//			return h.Coeff(i) == 0
//		})
//	if !hTest {
//		return false
//	}
//	// Reconstruct psi
//	At := A.Copy().AsMatrix().Transpose()
//	psi := math.NewMatrixFromDimensions(params.k, numSplits).PopulateRows(
//		func(mu int) math.Vector {
//			tmp := At.Copy().AsMatrix().MulVec(gamma.Row(mu))
//			return ens20.SplitInvNTT(tmp, numSplits, params.d, baseRing)
//		})
//	// Reconstruct the commitment to f
//	invk := math.NewModInt(int64(params.k), params.q).Inv()
//	tao := math.NewVectorFromSize(params.k).Populate(
//		func(mu int) math.RingElement {
//			tmp := math.NewVectorFromSize(params.k).Populate(
//				func(v int) math.RingElement {
//					tmp2 := math.NewVectorFromSize(numSplits).Populate(
//						func(j int) math.RingElement {
//							dec := math.NewOnePolynomial(baseRing).
//								Scale(u.Dot(gamma.Row(mu)).(*math.ModInt).Uint64()).
//								Neg()
//							return psi.Element(mu, j).Copy().
//								Mul(t.Element(j)).
//								Scale(uint64(params.d)).
//								Add(dec)
//						}).Sum()
//					return sig.Permute(int64(v), tmp2.(math.Polynomial))
//				}).Sum().(math.Polynomial)
//			return ens20.Lmu(mu, tmp, invk)
//		}).Sum()
//	// Verify the commitments
//	testResult := math.NewVectorFromSize(params.k).Populate(
//		func(i int) math.RingElement {
//			return math.NewVectorFromSize(params.k).Populate(
//				func(mu int) math.RingElement {
//					tmp := math.NewVectorFromSize(params.k).Populate(
//						func(v int) math.RingElement {
//							tmp2 := math.NewVectorFromSize(numSplits).Populate(
//								func(j int) math.RingElement {
//									index := (i - v) % params.k
//									return psi.Element(mu, j).Copy().
//										Scale(uint64(params.d)).
//										Mul(b.Row(j).Copy().AsVec().Dot(z.Row(index)))
//								}).Sum()
//							return sig.Permute(int64(v), tmp2.(math.Polynomial))
//						}).Sum().(math.Polynomial)
//					return ens20.Lmu(mu, tmp, invk)
//				}).Sum().Add(b.Row(numSplits + 1).Dot(z.Row(i)))
//		}).All(
//		func(lhs math.RingElement, i int) bool {
//			rhsAdd := sig.Permute(int64(i), c).Mul(
//				tao.Copy().Add(t.Element(numSplits + 1)).Copy().Add(h.Neg()))
//			rhs := vp.Element(i).Copy().Add(rhsAdd)
//			return lhs.Eq(rhs)
//		})
//	return testResult
//}
