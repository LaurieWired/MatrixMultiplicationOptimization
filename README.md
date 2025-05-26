# Matrix Multiplication Optimization

This repository demonstrates various matrix multiplication algorithms and their performance characteristics.

## Overview

This project includes:
1. **Naive vs. Strassen Algorithm Comparison** - Benchmarking traditional O(n³) multiplication against Strassen's O(n^2.807) algorithm
2. **Formula-based Matrix Multiplication** - Implementation of a parser that reads and applies optimized multiplication schemes from [mkauers/matrix-multiplication](https://github.com/mkauers/matrix-multiplication)
3. **Performance Analysis Tools** - Utilities to measure and compare the efficiency of different approaches

## Features

### 1. Algorithm Implementations

#### Naive Matrix Multiplication
- Standard triple-loop algorithm with O(n³) complexity
- Serves as a baseline for performance comparisons
- Used for verification of optimized algorithms

#### Strassen Algorithm (`strassen_multiply.cpp`)
- Divide-and-conquer approach reducing complexity to O(n^2.807)
- Implemented with configurable threshold for switching to naive multiplication
- Demonstrates practical speedup for larger matrices
- Automatic padding to power-of-2 dimensions for optimal recursive subdivision
- Threshold tuning capability (default: 64x64 matrices)

#### Formula-based Multiplication (`formula_multiply.cpp`)
- Parses `.m` formula files from the [mkauers/matrix-multiplication](https://github.com/mkauers/matrix-multiplication) repository
- Applies advanced multiplication schemes discovered through automated search techniques

### 2. Formula Parser

The `formula_multiply.cpp` program can parse and apply multiplication formulas in the format used by the research paper "Consequences of the Moosbauer-Poole Algorithms" (arXiv:2505.05896v1).

**Key capabilities:**
- Parses nested matrix structures from `.m` files
- Handles coefficient matrices for computing linear combinations
- Supports various matrix dimensions (4×5×6 through 6×6×7)
- Includes verification against naive multiplication

**Formula format:**
```
{{{{A_coeffs}}, {{B_coeffs}}, {{C_coeffs}}}, ...}
```
Each term represents a product in the formula: `C += (Σ A_coeffs ⊙ A) × (Σ B_coeffs ⊙ B) × C_coeffs`

### 3. Performance Results

Our tests show significant improvements using the formula-based approach:

| Matrix Size | Naive Ops | Formula Ops | Reduction % |
|-------------|-----------|-------------|-------------|
| 4×5×6       | 120       | 90          | 25.0%       |
| 4×5×7       | 140       | 104         | 25.7%       |
| 5×6×7       | 210       | 150         | 28.6%       |
| 6×6×7       | 252       | 183         | 27.4%       |

## Building the Project

### Prerequisites
- Formula files from [mkauers/matrix-multiplication](https://github.com/mkauers/matrix-multiplication) (included as submodule)

### Compilation

```bash
# Build all programs
make all

# Build specific targets
make formula_multiply
make strassen_multiply
make test_summary

# Clean build artifacts
make clean
```

## Usage

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

### Running Benchmarks

```bash
# Test all available matrix sizes
make test

# Display efficiency summary
make summary
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

## Project Structure

```
.
├── formula_multiply.cpp    # Formula parser and multiplier
├── strassen_multiply.cpp  # Strassen algorithm implementation
├── test_summary.cpp       # Efficiency analysis tool
├── test_all_sizes.sh      # Automated test script
├── Makefile              # Build configuration
├── README.md             # This file
```

## Technical Details

### Formula File Format

The `.m` files contain matrix multiplication schemes discovered through automated search. Each file represents a formula for multiplying matrices of specific dimensions, encoded as:
- Coefficient matrices for input matrix A (n×m)
- Coefficient matrices for input matrix B (m×p)
- Coefficient matrices for output matrix C (n×p)

### Implementation Notes

1. **Parser Design**: Uses recursive brace tracking to handle nested matrix structures
2. **Strassen Implementation**: Includes automatic padding and configurable threshold
3. **Optimization**: Skips zero coefficients to improve performance
4. **Verification**: Includes comparison with naive multiplication for correctness
5. **Precision**: Uses double-precision floating-point with configurable epsilon for comparisons

## Research Background

This implementation is based on the research paper:
- **"Consequences of the Moosbauer-Poole Algorithms"** by Manuel Kauers and Isaac Wood (arXiv:2505.05896v1)
- The formulas were discovered using flip graph search techniques
- Original formulas available at: [https://github.com/mkauers/matrix-multiplication](https://github.com/mkauers/matrix-multiplication)
