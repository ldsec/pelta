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
)

// DecomposeIntoTernary computes an integer vector v'S ternary (0, 1, 2) decomposition.
// Returns the ternary decomposition matrix and the basis vector.
func DecomposeIntoTernary(v rings.ZIntVector, logBeta int) (algebra.Matrix, rings.ZIntVector) {
	base := rings.NewZInt(3)
	basis := ComputeBasis(base, logBeta)
	decomposed := algebra.NewMatrixFromDimensions(v.Length(), logBeta)
	v.ForEach(
		func(el algebra.Element, i int) {
			// Get the basis representation of the value.
			basisRepr := IntoBasis(el.(rings.ZInt), base, logBeta).AsVec()
			decomposed.SetRow(i, basisRepr)
		})
	return decomposed, basis
}

// IntoBasis returns the given number'S representation under the given base, where the digits are
// encoded in the returned vector in an LSB manner.
func IntoBasis(num rings.ZInt, base rings.ZInt, logNum int) rings.ZIntVector {
	repr := make([]algebra.Element, 0, logNum)
	curr := num.Copy().(rings.ZInt)
	for i := 0; i < logNum; i++ {
		rem := curr.Copy().(rings.ZInt).Rem(base)
		curr.EuclideanDiv(base)
		repr = append(repr, rem)
	}
	return rings.NewZIntVec(algebra.NewVectorFromSlice(repr))
}

// ComputeBasis returns a vector of multiplicants S.t. (base^0, base^1, ..., base^(n-1)) that can be used for
// base decomposition of vectors.
func ComputeBasis(base rings.ZInt, n int) rings.ZIntVector {
	return rings.NewZIntVec(algebra.NewVectorFromSize(n).Populate(
		func(i int) algebra.Element {
			return base.Copy().Pow(uint64(i))
		}))
}

// GenerateNTTTransform computes and returns the integer NTT transformation matrix for the given base ring.
func GenerateNTTTransform(baseRing *ring.Ring, q *big.Int, logN uint64) algebra.Matrix {
	qInt := q.Uint64()
	w := ring.InvMForm(baseRing.NttPsi[0][baseRing.N>>1], qInt, baseRing.MredParams[0])
	mask := uint64(2*baseRing.N - 1)
	T := algebra.NewMatrixFromDimensions(baseRing.N, baseRing.N).PopulateRows(
		func(i int) algebra.Vector {
			twoirev := 2*utils.BitReverse64(uint64(i), logN) + 1
			return algebra.NewVectorFromSize(baseRing.N).Populate(
				func(j int) algebra.Element {
					gen := uint64(j) * twoirev & mask
					result := ring.ModExp(w, gen, qInt)
					return rings.NewZInt(int64(result))
				})
		})
	return T
}

// SaveNTTTransform generates the NTT transform matrix and saves it in a new file.
func SaveNTTTransform(baseRing *ring.Ring, q *big.Int, logN uint64) (algebra.Matrix, error) {
	fmt.Println("SaveNTTTransform: Generating the transform...")
	T := GenerateNTTTransform(baseRing, q, logN)
	fileName := fmt.Sprintf("NTT")
	file, err := os.Create(fileName)
	if err != nil {
		fmt.Println("couldn't generate the ntt transform file", err.Error())
		return algebra.Matrix{}, err
	}
	defer file.Close()
	fmt.Printf("SaveNTTTransform: Saving into file %s...\n", fileName)
	for i := 0; i < T.Length(); i++ {
		s := fmt.Sprintf("%d\n", T.MultiArray.Array[i].(rings.ZInt).Int64())
		file.WriteString(s)
	}
	fmt.Printf("SaveNTTTransform: Done\n")
	return T, nil
}

func LoadNTTTransform(baseRing *ring.Ring, q *big.Int) (algebra.Matrix, error) {
	fileName := fmt.Sprintf("NTT")
	file, err := os.Open(fileName)
	if err != nil {
		fmt.Println("couldn't open the ntt transformation file", err.Error())
		return algebra.Matrix{}, err
	}
	fmt.Printf("LoadNTTTransform: Loading the transform from file %s..\n.", fileName)
	rd := bufio.NewReader(file)
	T := algebra.NewMatrixFromDimensions(baseRing.N, baseRing.N)
	for i := 0; i < T.Length(); i++ {
		el, _ := rd.ReadString('\n')
		elI, _ := strconv.Atoi(el[:len(el)-1])
		T.SetElementAtIndex(i, rings.NewZqInt(int64(elI), q))
	}
	return T, nil
}
