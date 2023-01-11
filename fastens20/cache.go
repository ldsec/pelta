package fastens20

type Cache struct {
	values map[string]map[int64]uint64
}

func NewEmptyCache() *Cache {
	return &Cache{values: map[string]map[int64]uint64{}}
}

func (c *Cache) Get(name string, key int64, gen func() uint64) uint64 {
	_, ok := c.values[name]
	if !ok {
		c.values[name] = map[int64]uint64{}
	}
	_, ok = c.values[name][key]
	if !ok {
		c.values[name][key] = gen()
	}
	return c.values[name][key]
}
