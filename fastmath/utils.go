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
	return &loadedMatrix, nil
}
