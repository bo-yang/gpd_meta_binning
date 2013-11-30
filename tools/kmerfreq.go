// kmerfreq - count unique k-mers from a fasta file.

package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io"
	"os"
	"os/exec"
	"strings"
	"strconv"
)

// Readln returns a single line (without the ending \n)
// from the input buffered reader.
// An error is returned iff there is an error with the
// buffered reader.
func Readln(r *bufio.Reader) (string, error) {
	var (
		isPrefix bool  = true
		err      error = nil
		line, ln []byte
	)
	for isPrefix && err == nil {
		line, isPrefix, err = r.ReadLine()
		ln = append(ln, line...)
	}
	return string(ln), err
}

func main() {
	// Get options
	var input = flag.String("i", "", "input fasta file")
	var kmer = flag.Int("k", 16, "length of k-mer")
	flag.Parse()

	// Dynamically parse the input file name
	cmd := exec.Command("basename", *input)
	cmd.Stdin = strings.NewReader("some input")
	var name bytes.Buffer
	cmd.Stdout = &name
	err := cmd.Run()
	if err != nil {
		fmt.Println(err)
		return
	}
	var output = flag.String("o", fmt.Sprintf("%s.%d-mer", strings.TrimSpace(name.String()), *kmer), "output file")
	flag.Parse()

	// open input file
	fi, err := os.Open(*input)
	if err != nil {
		panic(err)
	}
	// close fi on exit and check for its returned error
	defer func() {
		if err := fi.Close(); err != nil {
			panic(err)
		}
	}()

	// open output file
	fo, err := os.Create(*output)
	if err != nil {
		panic(err)
	}
	// close fo on exit and check for its returned error
	defer func() {
		if err := fo.Close(); err != nil {
			panic(err)
		}
	}()

	// Read and count unique k-mers
	var read string = ""
	var km map[string]int  // k-mer counts
	km = make(map[string]int)

	r := bufio.NewReader(fi)
	for {
		line, err := Readln(r)
		if err == io.EOF {
			if !strings.HasPrefix(line, ">") {
				read += line
			}
			break // done
		} else if err != nil {
			panic(err) // error happens
		}
		if strings.HasPrefix(line, ">") {
			// Header of the read
			if len(read) != 0 {
				// count k-mers
				for i := 0; i < len(read)-*kmer+1; i++ {
					if km[read[i:i+*kmer]] != 0 {
						km[read[i:i+*kmer]]++
					} else {
						km[read[i:i+*kmer]] = 1
					}
				}
				read = ""
			}
		} else {
			// Read
			read += line
		}
	} // end of for

	// Don't forget the last read
	if len(read) != 0 {
		// count k-mers
		for i := 0; i < len(read)-*kmer+1; i++ {
			if km[read[i:i+*kmer]] != 0 {
				km[read[i:i+*kmer]]++
			} else {
				km[read[i:i+*kmer]] = 1
			}
		}
	}

	// Write results
	for k, v := range km {
		fo.WriteString(k + "\t" + strconv.Itoa(v) + "\n")
	}
}
