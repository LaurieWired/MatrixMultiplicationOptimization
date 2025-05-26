#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <iomanip>

// Define a type alias for matrices for convenience
using Matrix = std::vector<std::vector<double>>;

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

// Naive matrix multiplication
Matrix naiveMultiply(const Matrix& A, const Matrix& B) {
    if (A.empty() || B.empty() || A[0].size() != B.size()) {
        std::cerr << "Error: Incompatible matrix dimensions for naive multiplication." << std::endl;
        return {};
    }

    size_t n = A.size();    // Rows of A
    size_t m = B.size();    // Rows of B = Columns of A
    size_t p = B[0].size(); // Columns of B

    Matrix C = createMatrix(n, p);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < p; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < m; ++k) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
    return C;
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
     std::cout << "Naive Matrix Multiplication Benchmark" << std::endl;
    std::cout << "Matrix size: " << N << "x" << N << std::endl;


    Matrix A = createMatrix(N, N);
    Matrix B = createMatrix(N, N);

    fillMatrix(A);
    fillMatrix(B);

    // --- Time Naive Multiplication ---
    auto start_naive = std::chrono::high_resolution_clock::now();
    Matrix C_naive = naiveMultiply(A, B);
    auto end_naive = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> naive_duration = end_naive - start_naive;

    std::cout << "Naive multiplication time: " << naive_duration.count() << " ms" << std::endl;

    // if (N <= 8) {
    //     printMatrix(A, "A");
    //     printMatrix(B, "B");
    //     printMatrix(C_naive, "C (Naive)");
    // }
    if (C_naive.empty() && N > 0) {
        std::cerr << "Naive multiplication failed to produce a result." << std::endl;
    }


    return 0;
}