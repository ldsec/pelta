package crypto

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
)

type AjtaiCommitment struct {
	A     *fastmath.IntMatrix
	B     *fastmath.IntMatrix
	S     *fastmath.IntVec
	R     *fastmath.IntVec
	Kappa *fastmath.IntVec
	ComP  *fastmath.IntVec
	P     *big.Int
}

func NewAjtaiCommitment(s, r *fastmath.IntVec, comSize int, config Config) AjtaiCommitment {
	A := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	B := fastmath.NewRandomIntMatrix(comSize, s.Size(), config.P, config.BaseRing)
	// Compute kappa.
	comQ := A.MulVec(s).Add(B.MulVec(r))
	comP := comQ.Copy().Reduce(config.P)
	// comQ - comP = (As + Br) - [(As + Br) mod p]
	diff := comQ.Copy().Add(comP.Copy().Neg())
	// kappa := (comQ - comP) / p s.t. kappa * p is the difference between
	// pInv := big.NewInt(0).ModInverse(config.P, config.Q).Uint64()
	kappa := fastmath.NewIntVec(diff.Size(), config.BaseRing)
	kappa.Populate(func(i int) uint64 {
		return diff.Get(i) / config.P.Uint64()
	})
	return AjtaiCommitment{A, B, s, r, kappa, comP, config.P}
}

func (aj *AjtaiCommitment) EmbedIntoSIS(sis *SISProblem, config Config) {
	k := sis.S.Size()/config.D - 1
	l := aj.ComP.Size()
	// Extend u and s.
	uExtension := aj.ComP.Copy()
	sExtension := aj.R.Copy()
	sExtension.Append(aj.Kappa)
	// Pad the extensions.
	padLength := config.D - l
	padding := fastmath.NewIntVec(padLength, config.BaseRing)
	uExtension.Append(padding.Copy())
	sExtension.Append(padding)
	sis.U.Append(uExtension)
	sis.S.Append(sExtension)
	// Extend SIS A horizontally.
	zeroHorExt := fastmath.NewIntMatrix(config.D, 2*config.D, config.BaseRing)
	sis.A.ExtendCols(zeroHorExt)
	// Extend SIS A vertically.
	aExtensionParts := make([]fastmath.IntMatrix, k+3)
	aExtensionParts[0] = *aj.A.Copy()
	aExtensionParts[k+1] = *aj.B.Copy()
	// Pad the embedded A and B (# rows = l) to D rows
	zeroVerExt := fastmath.NewIntMatrix(config.D-l, config.D, config.BaseRing)
	aExtensionParts[0].ExtendRows(zeroVerExt.Copy())
	aExtensionParts[k+1].ExtendRows(zeroVerExt)
	for i := 0; i < k; i++ {
		aExtensionParts[i+1] = *fastmath.NewIntMatrix(config.D, config.D, config.BaseRing)
	}
	negP := config.Q.Uint64() - aj.P.Uint64()
	negPDiag := fastmath.NewIntMatrix(config.D, config.D, config.BaseRing)
	for i := 0; i < negPDiag.Rows(); i++ {
		negPDiag.Set(i, i, negP)
	}
	aExtensionParts[k+2] = *negPDiag
	aVertExtension := fastmath.NewIntMatrix(config.D, 7*config.D, config.BaseRing)
	for i := 0; i < aVertExtension.Rows(); i++ {
		aExtensionRow := aExtensionParts[0].RowView(i)
		for _, aExtPart := range aExtensionParts[1:] {
			aExtensionRow.Append(aExtPart.RowView(i))
		}
		aVertExtension.SetRow(i, *aExtensionRow)
	}
	sis.A.ExtendRows(aVertExtension)
}
