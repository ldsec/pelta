package rings

func min(a, b int64) int64 {
	if a <= b {
		return a
	}
	return b
}

func max(a, b int64) int64 {
	if a >= b {
		return a
	}
	return b
}

func Mod(a int, n int) int {
	a %= n
	if a < 0 {
		a += n
	}
	return a
}
