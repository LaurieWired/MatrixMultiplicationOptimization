#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <cmath>
#include <stdexcept>
#include <iomanip>

// Define a type alias for matrices for convenience
using Matrix = std::vector<std::vector<double>>;

// Threshold for switching from Strassen to naive (can be tuned)
const int STRASSEN_THRESHOLD = 64; // Power of 2 is often good

// Function to create an empty matrix
Matrix createMatrix(int rows, int cols) {
    return Matrix(rows, std::vector<double>(cols, 0.0));
}

// Function to fill a matrix with random values
void fillMatrix(Matrix& mat) {
    if (mat.empty() || mat[0].empty()) return;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distrib(1.0, 10.0);

    for (size_t i = 0; i < mat.size(); ++i) {
        for (size_t j = 0; j < mat[0].size(); ++j) {
            mat[i][j] = distrib(gen);
        }
    }
}

// Function to print a matrix (for debugging small matrices)
void printMatrix(const Matrix& mat, const std::string& name) {
    if (mat.empty()) {
        std::cout << name << " is empty." << std::endl;
        return;
    }
    std::cout << name << " (" << mat.size() << "x" << mat[0].size() << "):" << std::endl;
    for (const auto& row : mat) {
        for (double val : row) {
            std::cout << std::fixed << std::setprecision(2) << val << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// Naive matrix multiplication (used as base case for Strassen)
Matrix naiveMultiply(const Matrix& A, const Matrix& B) {
    if (A.empty() || B.empty() || A[0].size() != B.size()) {
        // This should ideally not happen if Strassen is called correctly
        // but good for standalone use or errors in recursion.
        std::cerr << "Error: Incompatible matrix dimensions for naive multiplication." << std::endl;
        return {};
    }
    size_t n = A.size();
    size_t m = B.size();
    size_t p = B[0].size();
    Matrix C = createMatrix(n, p);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < p; ++j) {
            for (size_t k = 0; k < m; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}


// Matrix addition
Matrix add(const Matrix& A, const Matrix& B) {
    size_t n = A.size();
    size_t m = A[0].size();
    Matrix C = createMatrix(n, m);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

// Matrix subtraction
Matrix subtract(const Matrix& A, const Matrix& B) {
    size_t n = A.size();
    size_t m = A[0].size();
    Matrix C = createMatrix(n, m);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

// Function to split a matrix into four submatrices
void splitMatrix(const Matrix& P, Matrix& P11, Matrix& P12, Matrix& P21, Matrix& P22) {
    size_t newSize = P.size() / 2;
    for (size_t i = 0; i < newSize; ++i) {
        for (size_t j = 0; j < newSize; ++j) {
            P11[i][j] = P[i][j];
            P12[i][j] = P[i][j + newSize];
            P21[i][j] = P[i + newSize][j];
            P22[i][j] = P[i + newSize][j + newSize];
        }
    }
}

// Function to join four submatrices into a single matrix
void joinMatrices(Matrix& P, const Matrix& P11, const Matrix& P12, const Matrix& P21, const Matrix& P22) {
    size_t newSize = P11.size();
    for (size_t i = 0; i < newSize; ++i) {
        for (size_t j = 0; j < newSize; ++j) {
            P[i][j] = P11[i][j];
            P[i][j + newSize] = P12[i][j];
            P[i + newSize][j] = P21[i][j];
            P[i + newSize][j + newSize] = P22[i][j];
        }
    }
}


// Strassen's algorithm recursive implementation
Matrix strassenRecursive(const Matrix& A, const Matrix& B) {
    size_t n = A.size();

    // Base case
    if (n <= STRASSEN_THRESHOLD) {
        return naiveMultiply(A, B);
    }

    // New size for submatrices
    size_t newSize = n / 2;

    // Create submatrices
    Matrix A11 = createMatrix(newSize, newSize);
    Matrix A12 = createMatrix(newSize, newSize);
    Matrix A21 = createMatrix(newSize, newSize);
    Matrix A22 = createMatrix(newSize, newSize);

    Matrix B11 = createMatrix(newSize, newSize);
    Matrix B12 = createMatrix(newSize, newSize);
    Matrix B21 = createMatrix(newSize, newSize);
    Matrix B22 = createMatrix(newSize, newSize);

    // Split matrices into submatrices
    splitMatrix(A, A11, A12, A21, A22);
    splitMatrix(B, B11, B12, B21, B22);

    // Compute the 7 products (M1 to M7 in some notations, P1 to P7 here)
    Matrix P1 = strassenRecursive(add(A11, A22), add(B11, B22));
    Matrix P2 = strassenRecursive(add(A21, A22), B11);
    Matrix P3 = strassenRecursive(A11, subtract(B12, B22));
    Matrix P4 = strassenRecursive(A22, subtract(B21, B11));
    Matrix P5 = strassenRecursive(add(A11, A12), B22);
    Matrix P6 = strassenRecursive(subtract(A21, A11), add(B11, B12));
    Matrix P7 = strassenRecursive(subtract(A12, A22), add(B21, B22));

    // Compute C's submatrices
    Matrix C11 = add(subtract(add(P1, P4), P5), P7);
    Matrix C12 = add(P3, P5);
    Matrix C21 = add(P2, P4);
    Matrix C22 = add(subtract(add(P1, P3), P2), P6); // P1 - P2 + P3 + P6

    // Join submatrices into the result matrix C
    Matrix C = createMatrix(n, n);
    joinMatrices(C, C11, C12, C21, C22);

    return C;
}

// Wrapper for Strassen: handles padding to power of 2
Matrix strassenMultiply(const Matrix& A_orig, const Matrix& B_orig) {
    if (A_orig.empty() || B_orig.empty() || A_orig[0].size() != B_orig.size()) {
        std::cerr << "Error: Incompatible matrix dimensions for Strassen multiplication." << std::endl;
        return {};
    }
    if (A_orig.size() != A_orig[0].size() || B_orig.size() != B_orig[0].size() || A_orig.size() != B_orig.size()){
        // Strassen as implemented here is simpler for square matrices
        // For rectangular, one would pad to a common square dimension, then extract.
        // Or use more complex block-Strassen variants.
        std::cerr << "Error: Strassen's algorithm as implemented here requires square matrices of the same size." << std::endl;
        std::cerr << "A: " << A_orig.size() << "x" << A_orig[0].size() << ", B: " << B_orig.size() << "x" << B_orig[0].size() << std::endl;
        // Fallback or error
        // return naiveMultiply(A_orig, B_orig);
        return {};
    }


    size_t n_orig = A_orig.size();
    if (n_orig == 0) return {};

    // Determine the size for padding (next power of 2)
    size_t m = 1;
    if (n_orig > 0 && (n_orig & (n_orig - 1)) != 0) { // if n_orig is not a power of 2
        m = static_cast<size_t>(std::pow(2, std::ceil(std::log2(static_cast<double>(n_orig)))));
    } else {
        m = n_orig; // Already a power of 2 (or 1)
    }
    if (m == 0 && n_orig > 0) m = n_orig; // Handle case where n_orig is 1, log2(1)=0, pow(2,0)=1
                                         // or if n_orig is tiny causing m to be 0.
    if (m < STRASSEN_THRESHOLD && (m & (m-1)) !=0 ) { // if padded size is small and not power of 2
        // just pad to threshold if it's a power of 2, or next power of 2
         size_t next_pow2_threshold = STRASSEN_THRESHOLD;
         if ((STRASSEN_THRESHOLD & (STRASSEN_THRESHOLD-1)) != 0) { // if threshold not power of 2
            next_pow2_threshold = static_cast<size_t>(std::pow(2, std::ceil(std::log2(static_cast<double>(STRASSEN_THRESHOLD)))));
         }
         m = std::max(m, next_pow2_threshold);
         if (n_orig > m) m = n_orig; // Ensure we don't shrink
    }


    Matrix A_padded = createMatrix(m, m);
    Matrix B_padded = createMatrix(m, m);

    for (size_t i = 0; i < n_orig; ++i) {
        for (size_t j = 0; j < n_orig; ++j) {
            A_padded[i][j] = A_orig[i][j];
            B_padded[i][j] = B_orig[i][j];
        }
    }

    Matrix C_padded = strassenRecursive(A_padded, B_padded);

    // Extract the original size result from the padded result
    Matrix C_final = createMatrix(n_orig, n_orig);
    for (size_t i = 0; i < n_orig; ++i) {
        for (size_t j = 0; j < n_orig; ++j) {
            C_final[i][j] = C_padded[i][j];
        }
    }

    return C_final;
}


int main(int argc, char* argv[]) {
    int N = 512; // Default matrix size (N x N)
    if (argc > 1) {
        try {
            N = std::stoi(argv[1]);
             if (N <= 0) {
                std::cerr << "Matrix size must be positive. Using default " << 512 << "." << std::endl;
                N = 512;
            }
        } catch (const std::exception& e) {
            std::cerr << "Invalid argument for N. Using default " << 512 << ". Error: " << e.what() << std::endl;
            N = 512;
        }
    }
    std::cout << "Strassen's Matrix Multiplication Benchmark" << std::endl;
    std::cout << "Matrix size: " << N << "x" << N << std::endl;
    std::cout << "Strassen threshold: " << STRASSEN_THRESHOLD << std::endl;

    Matrix A = createMatrix(N, N);
    Matrix B = createMatrix(N, N);

    fillMatrix(A);
    fillMatrix(B);

    // --- Time Strassen Multiplication ---
    auto start_strassen = std::chrono::high_resolution_clock::now();
    Matrix C_strassen = strassenMultiply(A, B);
    auto end_strassen = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> strassen_duration = end_strassen - start_strassen;

    std::cout << "Strassen multiplication time: " << strassen_duration.count() << " ms" << std::endl;

    // if (N <= 8) {
    //     printMatrix(A, "A");
    //     printMatrix(B, "B");
    //     printMatrix(C_strassen, "C (Strassen)");
    // }
    if (C_strassen.empty() && N > 0) {
        std::cerr << "Strassen multiplication failed to produce a result." << std::endl;
    }

    return 0;
}