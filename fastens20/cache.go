package fastens20

import (
	"github.com/ldsec/codeBase/commitment/fastmath"
	"github.com/ldsec/codeBase/commitment/logging"
)

type Cache struct {
	SigmaExpCache map[int64]uint64
}

func NewEmptyCache() *Cache {
	return &Cache{
		SigmaExpCache: map[int64]uint64{},
	}
}

func (c *Cache) Build(k int, sig fastmath.Automorphism) {
	e := logging.LogExecStart("Cache", "building")
	defer e.LogExecEnd()
	for i := 0; i < k; i++ {
		c.SigmaExpCache[int64(i)] = sig.Exponent(int64(i))
		c.SigmaExpCache[int64(-i)] = sig.Exponent(int64(-i))
	}
}
