package rings

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/math/algebra"
)

// ZInt represents an integer in the ring Z, optionally in Zq
type ZInt struct {
	Value   *big.Int
	Modulus *big.Int
}

func NewZeroZInt() ZInt {
	newValue := big.NewInt(0)
	i := ZInt{newValue, nil}
	return i
}

func NewOneZInt() ZInt {
	newValue := big.NewInt(0)
	i := ZInt{newValue, nil}
	return i
}

func NewOneZqInt(mod *big.Int) ZInt {
	newValue := big.NewInt(1)
	i := ZInt{newValue, mod}
	return i
}

func NewZInt(value int64) ZInt {
	i := ZInt{big.NewInt(value), nil}
	return i
}

func NewZqInt(value int64, mod *big.Int) ZInt {
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
	return NewZqInt(m.Value.Int64(), m.Modulus)
}

func (m ZInt) Pow(exp uint64) algebra.Element {
	m.Value = m.Value.Exp(m.Value, big.NewInt(int64(exp)), m.Modulus)
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

func (m ZInt) Rem(el ZInt) ZInt {
	m.Value.Rem(m.Value, el.Value)
	return m
}

func (m ZInt) EuclideanDiv(el ZInt) ZInt {
	m.Value.Div(m.Value, el.Value)
	return m
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
