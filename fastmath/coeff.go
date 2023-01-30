package fastmath

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// Coeff represents a coefficient of a polynomial over all levels.
type Coeff []uint64

func NewZeroCoeff(length int) Coeff {
	return make([]uint64, length)
}

func NewCoeffFromUint64(v uint64, mods []uint64) Coeff {
	c := make([]uint64, len(mods))
	for lvl, mod := range mods {
		c[lvl] = v % mod
	}
	return c
}

func NewCoeffFromBigInt(v *big.Int, mods []uint64) Coeff {
	c := make([]uint64, len(mods))
	for lvl, mod := range mods {
		c[lvl] = big.NewInt(0).Mod(v, big.NewInt(int64(mod))).Uint64()
	}
	return c
}

// IsZero returns true iff this coefficient is zero.
func (c Coeff) IsZero() bool {
	return c[0] == 0
}

// AsBigInt returns the big int representation of this coefficient.
func (c Coeff) AsBigInt(q *big.Int) *big.Int {
	acc := big.NewInt(0)
	for _, cl := range c {
		acc = acc.Mul(acc, big.NewInt(int64(cl)))
	}
	acc.Mod(acc, q)
	return acc
}

// Neg negates the coefficient across levels.
func (c Coeff) Neg(mods []uint64) Coeff {
	for i := range c {
		c[i] = mods[i] - c[i]
	}
	return c
}

func (c Coeff) Add(c2 Coeff, mods []uint64) Coeff {
	for i := range c {
		c[i] = (c[i] + c2[i]) % mods[i]
	}
	return c
}

// ExtendAsPoly returns a new polynomial with all of its coefficients set to this value.
func (c Coeff) ExtendAsPoly(baseRing *ring.Ring) *Poly {
	p := NewPoly(baseRing)
	for i := 0; i < p.N(); i++ {
		p.SetCoeff(i, c)
	}
	return p
}
