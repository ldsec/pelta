package math

// Utility functions that should've been in the std in the first place. Thank you Rob Pike.

func (m *MultiArray) min(a, b int) int {
	if a <= b {
		return a
	}
	return b
}

func (m *MultiArray) max(a, b int) int {
	if a >= b {
		return a
	}
	return b
}

func (m *MultiArray) abs(a int) int {
	if a < 0 {
		return -a
	}
	return a
}
