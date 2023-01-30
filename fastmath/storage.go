package fastmath

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/ldsec/codeBase/commitment/logging"
	"github.com/tuneinsight/lattigo/v4/ring"
)

// PersistentIntMatrix either loads the matrix in the file with `name` or generates and saves it using the given generator.
func PersistentIntMatrix(name string, generator func() *IntMatrix, baseRing *ring.Ring) *IntMatrix {
	procName := fmt.Sprintf("PersistentIntMatrix(%s)", name)
	var M *IntMatrix
	if _, err := os.Stat(name); errors.Is(err, os.ErrNotExist) {
		M = logging.LogExecution(procName, "generation", func() interface{} { return generator() }).(*IntMatrix)
		logging.LogShortExecution(procName, "saving", func() interface{} {
			err := SaveIntMatrix(M, name)
			if err != nil {
				logging.Log(procName, fmt.Sprintf("couldn't save matrix %s", err.Error()))
			}
			return nil
		})
	} else {
		M = logging.LogShortExecution(procName, "loading", func() interface{} {
			M, err := LoadIntMatrix(name, baseRing)
			if err != nil {
				logging.Log(procName, fmt.Sprintf("couldn't load matrix %s", err.Error()))
				M = logging.LogExecution(procName, "generation", func() interface{} { return generator() }).(*IntMatrix)
			}
			return M
		}).(*IntMatrix)
	}
	return M
}

// PersistentIntVec either loads the vector in the file with `name` or generates and saves it using the given generator.
func PersistentIntVec(name string, generator func() *IntVec, baseRing *ring.Ring) *IntVec {
	procName := fmt.Sprintf("PersistentIntVec(%s)", name)
	var M *IntVec
	if _, err := os.Stat(name); errors.Is(err, os.ErrNotExist) {
		M = logging.LogExecution(procName, "generation", func() interface{} { return generator() }).(*IntVec)
		logging.LogShortExecution(procName, "saving", func() interface{} {
			err := SaveIntVec(M, name)
			if err != nil {
				logging.Log(procName, fmt.Sprintf("couldn't save vector %s", err.Error()))
			}
			return nil
		})
	} else {
		M = logging.LogShortExecution(procName, "loading", func() interface{} {
			M, err := LoadIntVec(name, baseRing)
			if err != nil {
				logging.Log(procName, fmt.Sprintf("couldn't load matrix %s", err.Error()))
				M = logging.LogExecution(procName, "generation", func() interface{} { return generator() }).(*IntVec)
			}
			return M
		}).(*IntVec)
	}
	return M
}

// SaveIntVec saves the given vector into the specified file.
func SaveIntVec(t *IntVec, name string) error {
	m := NewIntMatrix(1, t.Size(), t.baseRing)
	m.SetRow(0, t)
	return SaveIntMatrix(m, name)
}

// LoadIntVec loads the vector saved in the given file.
func LoadIntVec(name string, baseRing *ring.Ring) (*IntVec, error) {
	m, err := LoadIntMatrix(name, baseRing)
	if err != nil {
		return nil, err
	}
	return m.RowView(0), nil
}

// SaveIntMatrix saves the given matrix into the specified file.
func SaveIntMatrix(t *IntMatrix, name string) error {
	file, err := os.Create(name)
	if err != nil {
		return err
	}
	defer file.Close()
	rowStrings := make([]string, 0, t.Rows())
	for _, row := range t.rows {
		colStrings := make([]string, 0, row.Size())
		for col := 0; col < row.Size(); col++ {
			colStrings = append(colStrings, unparseCoeff(row.GetCoeff(col)))
		}
		// Separate column values by a comma.
		rowStr := strings.Join(colStrings, ",")
		rowStrings = append(rowStrings, rowStr)
	}
	// Separate the rows by a newline.
	s := strings.Join(rowStrings, "\n")
	_, err = file.WriteString(s)
	return err
}

// LoadIntMatrix loads the matrix saved in the given file.
func LoadIntMatrix(name string, baseRing *ring.Ring) (*IntMatrix, error) {
	file, err := os.Open(name)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	loadedRows := make([][]Coeff, 0)
	bufReader := bufio.NewReader(file)
	for {
		line, err := bufReader.ReadString('\n')
		// Handle unexpected error.
		if err != nil && !errors.Is(err, io.EOF) {
			return nil, err
		}
		// Discard the newline character.
		if line[len(line)-1] == '\n' {
			line = line[:len(line)-1]
		}
		// Read in the column values
		colStrs := strings.Split(line, ",")
		matrixRow := make([]Coeff, 0, len(colStrs))
		for _, colStr := range colStrs {
			col, err := parseCoeff(colStr)
			if err != nil {
				return nil, err
			}
			matrixRow = append(matrixRow, col)
		}
		loadedRows = append(loadedRows, matrixRow)
		// Break on the last line.
		if err != nil && errors.Is(err, io.EOF) {
			break
		}
	}
	loadedMatrix := NewIntMatrixFromCoeffSlice(loadedRows, baseRing)
	return loadedMatrix, nil
}

func parseCoeff(s string) (Coeff, error) {
	vals := make([]uint64, 0)
	for _, vstr := range strings.Split(s, ";") {
		v, err := strconv.ParseUint(vstr, 10, 64)
		if err != nil {
			return nil, err
		}
		vals = append(vals, v)
	}
	return vals, nil
}

func unparseCoeff(c Coeff) string {
	valStrs := make([]string, len(c))
	for i, cl := range c {
		valStrs[i] = fmt.Sprintf("%d", cl)
	}
	return fmt.Sprintf("%s", strings.Join(valStrs, ";"))
}
