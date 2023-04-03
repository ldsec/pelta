package fastens20

import (
	"math/big"

	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
)

// Cache is a cache that contains the precomputed values to be used
// in the ENS20 protocol execution.
type Cache struct {
	SigmaExpCache map[int64]uint64
	InvK          uint64
}

// NewEmptyCache creates an empty cache.
func NewEmptyCache() *Cache {
	return &Cache{
		SigmaExpCache: map[int64]uint64{},
	}
}

// Build builds the protocol cache.
func (c *Cache) Build(k int, sig fastmath.Automorphism, q *big.Int) {
	e := logging.LogExecStart("Cache", "building")
	defer e.LogExecEnd()
	for i := 0; i < k; i++ {
		c.SigmaExpCache[int64(i)] = sig.Exponent(int64(i))
		c.SigmaExpCache[int64(-i)] = sig.Exponent(int64(-i))
	}
	c.InvK = big.NewInt(0).ModInverse(big.NewInt(int64(k)), q).Uint64()
}
