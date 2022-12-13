package fastmath

func AcceptIntSample(z *IntVec, samplingBound uint64) bool {
	return z.Max() < samplingBound
}

func AcceptPolyVecSample(z *PolyVec, samplingBound uint64) bool {
	return z.All(func(_ int, p *Poly) bool {
		return p.Max(0) < samplingBound
	})
}
