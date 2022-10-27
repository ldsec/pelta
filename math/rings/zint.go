package rings

import (
	"github.com/ldsec/codeBase/commitment/math/algebra"
	"math/big"
)

// ZInt represents an integer in the ring Z, optionally in Zq
type ZInt struct {
	Value   *big.Int
	Modulus *big.Int
}

func NewZeroInt() ZInt {
	newValue := big.NewInt(0)
	i := ZInt{newValue, nil}
	return i
}

func NewOneInt() ZInt {
	newValue := big.NewInt(0)
	i := ZInt{newValue, nil}
	return i
}

func NewOneModInt(mod *big.Int) ZInt {
	newValue := big.NewInt(1)
	i := ZInt{newValue, mod}
	return i
}

func NewInt(value int64) ZInt {
	i := ZInt{big.NewInt(value), nil}
	return i
}

func NewModInt(value int64, mod *big.Int) ZInt {
	i := ZInt{big.NewInt(value), mod}
	i.Reduce()
	return i
}

func (m ZInt) Add(q algebra.Element) algebra.Element {
	m.Value.Add(m.Value, q.(ZInt).Value)
	m.Reduce()
	return m
}

func (m ZInt) Sub(q algebra.Element) algebra.Element {
	return m.Add(q.Copy().Neg())
}

func (m ZInt) Mul(q algebra.Element) algebra.Element {
	m.Value.Mul(m.Value, q.(ZInt).Value)
	m.Reduce()
	return m
}

func (m ZInt) MulAdd(q algebra.Element, out algebra.Element) {
	// Keep m * q
	tmp := m.Copy()
	tmp.Mul(q)
	out.Add(tmp)
}

func (m ZInt) Neg() algebra.Element {
	m.Value.Neg(m.Value)
	m.Reduce()
	return m
}

func (m ZInt) Zero() algebra.Element {
	m.Value.Set(big.NewInt(0))
	return m
}

func (m ZInt) One() algebra.Element {
	m.Value.Set(big.NewInt(1))
	return m
}

func (m ZInt) Copy() algebra.Element {
	return NewModInt(m.Value.Int64(), m.Modulus)
}

func (m ZInt) Pow(exp uint64) algebra.Element {
	out := m.Copy().One().(ZInt)
	for i := uint64(0); i < exp; i++ {
		out.Mul(m)
	}
	m.Value = out.Value
	return m
}

func (m ZInt) Scale(factor uint64) algebra.Element {
	m.Value.Mul(m.Value, big.NewInt(int64(factor)))
	m.Reduce()
	return m
}

func (m ZInt) Eq(el algebra.Element) bool {
	m.Reduce()
	el.(ZInt).Reduce()
	return el.(ZInt).Value.Cmp(m.Value) == 0
}

func (m ZInt) String() string {
	m.Reduce()
	return "Int{" + m.Value.String() + "}"
}

// Inv sets this integer to its multiplicative inverse and returns the result.
func (m ZInt) Inv() ZInt {
	if m.Modulus != nil {
		m.Value.ModInverse(m.Value, m.Modulus)
	}
	return m
}

// Reduce reduces this integer to its mod ring.
func (m ZInt) Reduce() ZInt {
	if m.Modulus != nil {
		m.Value.Mod(m.Value, m.Modulus)
	}
	return m
}

// Int64 returns the int64 representation of this integer.
func (m ZInt) Int64() int64 {
	m.Reduce()
	return m.Value.Int64()
}

// Uint64 returns the int64 representation of this integer.
func (m ZInt) Uint64() uint64 {
	m.Reduce()
	if m.Modulus == nil {
		panic("cannot convert Z int to uint!!")
	}
	return m.Value.Uint64()
}

// Uint64WithMod returns the int64 representation of this integer.
func (m ZInt) Uint64WithMod(mod *big.Int) uint64 {
	m.Reduce()
	repr := big.NewInt(m.Value.Int64())
	return repr.Mod(repr, mod).Uint64()
}
