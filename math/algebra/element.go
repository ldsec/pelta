package algebra

// Element represents an element of a ring.
type Element interface {
	Add(q Element) Element         // p = p + q
	Sub(q Element) Element         // p = p - q
	Mul(q Element) Element         // p = p * q
	MulAdd(q Element, out Element) // out += p * q
	Neg() Element                  // p = -p
	Pow(exp uint64) Element        // p = p^exp
	Scale(factor uint64) Element   // p = factor*p
	Zero() Element                 // Converts into additive identity
	One() Element                  // Converts into multiplicative identity
	Copy() Element                 // Returns a copy of the element
	Eq(Element) bool               // Returns true if the two elements are equal
	String() string
}
