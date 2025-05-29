# Matrix Multiplication Optimization

This repository demonstrates algorithms for matrix multiplication to optimize both performance and rank.

## Overview

This project includes:
1. **Naive vs. Strassen Algorithm Comparison** - Benchmarking traditional O(n³) multiplication against Strassen's O(n^2.807) algorithm
2. **Formula-based Matrix Multiplication** - Parser that reads and applies optimized multiplication schemes from [mkauers/matrix-multiplication](https://github.com/mkauers/matrix-multiplication)

## Usage

### Prerequisites
- Formula files from [mkauers/matrix-multiplication](https://github.com/mkauers/matrix-multiplication)

### Running Strassen vs. Naive Benchmark

```bash
# Run benchmark with default 512x512 matrices
./strassen_multiply

# Run with custom size (must be square matrices)
./strassen_multiply 1024
```

### Running the Formula-based Multiplication

```bash
# Basic usage
./formula_multiply <formula_file.m> <n> <m> <p>

# Example: 4×5×6 matrix multiplication
./formula_multiply matrix-multiplication/456/k0e7ba384e845ae2.m 4 5 6

# With verification against naive algorithm
./formula_multiply matrix-multiplication/456/k0e7ba384e845ae2.m 4 5 6 --verify
```

### Available Make Targets

- `make all` - Build all executables
- `make formula_multiply` - Build the formula parser/multiplier
- `make strassen_multiply` - Build Strassen algorithm benchmark
- `make test_summary` - Build the efficiency analysis tool
- `make test` - Run tests on all matrix sizes
- `make summary` - Display efficiency comparison table
- `make clean` - Remove build artifacts
- `make help` - Show all available targets
