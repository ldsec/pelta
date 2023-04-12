package fastmath

import (
	"math/big"
	"sync"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// Coeff represents a coefficient of a polynomial over all levels.
type Coeff []uint64

// NewZeroCoeff creates a zero coefficient.
func NewZeroCoeff(length int) Coeff {
	return make([]uint64, length)
}

// NewCoeffFromUint64 creates a new coefficient from the given uint64 by taking its mod for every level mod
// given in `mods`.
func NewCoeffFromUint64(v uint64, mods []uint64) Coeff {
	c := make([]uint64, len(mods))
	var wg sync.WaitGroup
	wg.Add(len(mods))
	for lvl, mod := range mods {
		func(lvl int, mod uint64) {
			c[lvl] = v % mod
			wg.Done()
		}(lvl, mod)
	}
	wg.Wait()
	return c
}

// NewCoeffFromBigInt creates a new coefficient from the given bigint by taking its mod for every level mod
// given in `mods`.
func NewCoeffFromBigInt(v *big.Int, mods []uint64) Coeff {
	c := make([]uint64, len(mods))
	var wg sync.WaitGroup
	wg.Add(len(mods))
	for lvl, mod := range mods {
		func(lvl int, mod uint64) {
			qlvl := big.NewInt(int64(mod))
			c[lvl] = big.NewInt(0).Mod(v, qlvl).Uint64()
			wg.Done()
		}(lvl, mod)
	}
	wg.Wait()
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
	var wg sync.WaitGroup
	wg.Add(len(c))
	for i := range c {
		func(i int) {
			c[i] = mods[i] - c[i]
			wg.Done()
		}(i)
	}
	wg.Wait()
	return c
}

// Add adds two coefficients.
func (c Coeff) Add(c2 Coeff, mods []uint64) Coeff {
	var wg sync.WaitGroup
	wg.Add(len(c))
	for i := range c {
		func(i int) {
			c[i] = (c[i] + c2[i]) % mods[i]
			wg.Done()
		}(i)
	}
	wg.Wait()
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
