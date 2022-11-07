package rlwe

import (
	"bufio"
	"fmt"
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"github.com/ldsec/codeBase/commitment/math/rings"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
	"math/big"
	"os"
	"strconv"
	"strings"
)

// NTTTransformMatrix encodes the NTT transformation, represented as a vector of polynomials.
type NTTTransformMatrix struct {
	rings.PolyVector
	logN int
}

// Extended returns an extension of the transformation which includes the multiplication with
// a polynomial at the NTT domain.
func (m NTTTransformMatrix) Extended(poly rings.Polynomial) NTTTransformMatrix {
	copiedTransform := m.Copy().AsVec()
	Tpoly := poly.Copy().(rings.Polynomial).NTT()
	// Scale each element by the corresponding coefficient.
	copiedTransform.Map(
		func(el algebra.Element, i int) algebra.Element {
			return el.Scale(Tpoly.Coeff(i))
		})
	// Return the scaled transform.
	return NTTTransformMatrix{rings.NewPolyVec(copiedTransform), m.logN}
}

func (m NTTTransformMatrix) Transposed() NTTTransformMatrix {
	transposed := make([]algebra.Element, m.Length())
	for i := 0; i < m.Length(); i++ {
		poly := rings.NewZeroPolynomial(m.Array[0].(rings.Polynomial).BaseRing)
		for j, tj := range m.Array {
			poly.SetCoefficient(j, tj.(rings.Polynomial).Coeff(i))
		}
		transposed[i] = poly
	}
	return NTTTransformMatrix{rings.NewPolyVec(algebra.NewVectorFromSlice(transposed)), m.logN}
}

// Scaled returns an extension of the transformation which includes the multiplication with a scalar.
func (m NTTTransformMatrix) Scaled(factor uint64) NTTTransformMatrix {
	copiedTransform := m.Copy().Scale(factor).AsVec()
	// Return the scaled transform.
	return NTTTransformMatrix{rings.NewPolyVec(copiedTransform), m.logN}
}

func (m NTTTransformMatrix) Apply(p rings.Polynomial) rings.Polynomial {
	outPoly := rings.NewZeroPolynomial(p.BaseRing)
	for i := 0; i < p.BaseRing.N; i++ {
		// Get a temporary vector to store T[i]s
		tmp := rings.NewZeroPolynomial(p.BaseRing)
		// Mult s in poly form by T[i] coeff-wise
		p.BaseRing.MulCoeffs(p.Ref, m.Array[i].(rings.Polynomial).Ref, tmp.Ref)
		// Take the sum of all coefficients
		for j := 0; j < m.logN; j++ {
			tmp2 := rings.NewZeroPolynomial(p.BaseRing)
			p.BaseRing.Shift(tmp.Ref, 1<<j, tmp2.Ref)
			p.BaseRing.Add(tmp.Ref, tmp2.Ref, tmp.Ref)
		}
		// Copy the 1-st slot of tmp into the i-th slot of u.
		outPoly.SetCoefficient(i, tmp.Coeff(0))
	}
	outPoly.Ref.IsNTT = true
	return outPoly
}

// GenerateNTTTransform computes and returns the integer NTT transformation matrix for the given base ring.
func GenerateNTTTransform(baseRing *ring.Ring, q *big.Int, logN int) NTTTransformMatrix {
	qInt := q.Uint64()
	w := ring.InvMForm(baseRing.NttPsi[0][baseRing.N>>1], qInt, baseRing.MredParams[0])
	mask := uint64(2*baseRing.N - 1)
	T := algebra.NewVectorFromSize(baseRing.N).Populate(
		func(i int) algebra.Element {
			twoirev := 2*utils.BitReverse64(uint64(i), uint64(logN)) + 1
			poly := rings.NewZeroPolynomial(baseRing)
			for j := 0; j < baseRing.N; j++ {
				gen := uint64(j) * twoirev & mask
				poly.SetCoefficient(j, ring.ModExp(w, gen, qInt))
			}
			return poly
		})
	return NTTTransformMatrix{rings.NewPolyVec(T), logN}
}

// SaveNTTTransform generates the NTT transform matrix and saves it in a new file.
func SaveNTTTransform(baseRing *ring.Ring, q *big.Int, logN int) (NTTTransformMatrix, error) {
	fmt.Println("SaveNTTTransform: Generating the transform...")
	T := GenerateNTTTransform(baseRing, q, logN)
	fileName := fmt.Sprintf("NTT")
	file, err := os.Create(fileName)
	if err != nil {
		fmt.Println("couldn't generate the ntt transform file", err.Error())
		return NTTTransformMatrix{}, err
	}
	defer file.Close()
	fmt.Printf("SaveNTTTransform: Saving into file %s...\n", fileName)
	for i := 0; i < T.Length(); i++ {
		for j := 0; j < baseRing.N; j++ {
			s := fmt.Sprintf("%d,", int64(T.Element(i).(rings.Polynomial).Coeff(j)))
			file.WriteString(s)
		}
		file.WriteString("\n")
	}
	fmt.Printf("SaveNTTTransform: Done\n")
	return T, nil
}

func LoadNTTTransform(baseRing *ring.Ring, logN int) (NTTTransformMatrix, error) {
	fileName := fmt.Sprintf("NTT")
	file, err := os.Open(fileName)
	if err != nil {
		fmt.Println("couldn't open the ntt transformation file", err.Error())
		return NTTTransformMatrix{}, err
	}
	fmt.Printf("LoadNTTTransform: Loading the transform from file %s...\n", fileName)
	rd := bufio.NewReader(file)
	T := algebra.NewVectorFromSize(baseRing.N)
	for i := 0; i < T.Length(); i++ {
		// Read the polynomial.
		polyString, _ := rd.ReadString('\n')
		coeffStrings := strings.Split(polyString[:len(polyString)-1], ",")
		poly := rings.NewZeroPolynomial(baseRing)
		for j, coeffString := range coeffStrings[:len(coeffStrings)-1] {
			coeff, _ := strconv.Atoi(coeffString)
			poly.SetCoefficient(j, uint64(coeff))
		}
		// Save the polynomial.
		T.SetElementAtIndex(i, poly)
	}
	fmt.Printf("LoadNTTTransform: Done.\n")
	return NTTTransformMatrix{rings.NewPolyVec(T), logN}, nil
}
