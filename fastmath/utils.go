package fastmath

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// PersistentIntMatrix either loads the matrix in the file with `name` or generates and saves it using the generator.
func PersistentIntMatrix(name string, generator func() *IntMatrix, baseRing *ring.Ring) *IntMatrix {
	var M *IntMatrix
	if _, err := os.Stat(name); errors.Is(err, os.ErrNotExist) {
		fmt.Printf("PersistentIntMatrix: couldn't find the file %s, generating...", name)
		M = generator()
		fmt.Printf("done\n")
		fmt.Printf("PersistentIntMatrix: saving...")
		err := SaveIntMatrix(M, name)
		fmt.Printf("done\n")
		if err != nil {
			fmt.Printf("PersistentIntMatrix: couldn't save matrix %s: %s\n", name, err.Error())
		}
	} else {
		fmt.Printf("PersistentIntMatrix: loading %s...", name)
		M, err = LoadIntMatrix(name, baseRing)
		fmt.Printf("done\n")
		if err != nil {
			fmt.Printf("PersistentIntMatrix: couldn't load matrix %s, regenerating: %s", name, err.Error())
			M = generator()
			fmt.Printf("done\n")
		}
	}
	return M
}

// PersistentIntVec either loads the vector in the file with `name` or generates and saves it using the generator.
func PersistentIntVec(name string, generator func() *IntVec, baseRing *ring.Ring) *IntVec {
	var v *IntVec
	if _, err := os.Stat(name); errors.Is(err, os.ErrNotExist) {
		fmt.Printf("PersistentIntVec: couldn't find the file %s, generating...", name)
		v = generator()
		fmt.Printf("done\n")
		fmt.Printf("PersistentIntVec: saving...")
		err := SaveIntVec(v, name)
		fmt.Printf("done\n")
		if err != nil {
			fmt.Printf("PersistentIntVec: couldn't save matrix %s: %s\n", name, err.Error())
		}
	} else {
		fmt.Printf("PersistentIntVec: loading %s...", name)
		v, err = LoadIntVec(name, baseRing)
		fmt.Printf("done\n")
		if err != nil {
			fmt.Printf("PersistentIntVec: couldn't load matrix %s, regenerating: %s", name, err.Error())
			v = generator()
			fmt.Printf("done\n")
		}
	}
	return v
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
			colStr := fmt.Sprintf("%d", row.Get(col))
			colStrings = append(colStrings, colStr)
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
	loadedRows := make([][]uint64, 0)
	bufReader := bufio.NewReader(file)
	for {
		line, err := bufReader.ReadString('\n')
		// Handle unexpected error.
		if err != nil && !errors.Is(err, io.EOF) {
			fmt.Printf("encountered an unexpected error while loading the int matrix: %s\n", err.Error())
			return nil, err
		}
		// Discard the newline character.
		if line[len(line)-1] == '\n' {
			line = line[:len(line)-1]
		}
		// Read in the column values
		colStrs := strings.Split(line, ",")
		matrixRow := make([]uint64, 0, len(colStrs))
		for _, colStr := range colStrs {
			col, err := strconv.ParseUint(colStr, 10, 64)
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
	loadedMatrix := NewIntMatrixFromSlice(loadedRows, baseRing)
	return loadedMatrix, nil
}
