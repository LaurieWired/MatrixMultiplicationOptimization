#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <random>
#include <chrono>

using Matrix = std::vector<std::vector<double>>;

// --- Matrix Utilities ---
Matrix createMatrix(int rows, int cols) {
    if (rows < 0 || cols < 0) {
        throw std::runtime_error("Matrix dimensions cannot be negative.");
    }
    return Matrix(static_cast<size_t>(rows), std::vector<double>(static_cast<size_t>(cols), 0.0));
}

void fillMatrixRandom(Matrix& mat) {
    if (mat.empty() || mat[0].empty()) return;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distrib(1.0, 5.0);
    for (size_t i = 0; i < mat.size(); ++i) {
        for (size_t j = 0; j < mat[0].size(); ++j) {
            mat[i][j] = distrib(gen);
        }
    }
}

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

// Naive matrix multiplication for verification
Matrix naiveMultiply(const Matrix& A, const Matrix& B) {
    if (A.empty() || B.empty() || A[0].size() != B.size()) {
        throw std::runtime_error("Error: Incompatible matrix dimensions for naive multiplication.");
    }

    size_t n = A.size();    // Rows of A
    size_t m_common = B.size();    // Rows of B = Columns of A
    size_t p = B[0].size(); // Columns of B

    Matrix C = createMatrix(n, p);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < p; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < m_common; ++k) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
    return C;
}

bool compareMatrices(const Matrix& M1, const Matrix& M2, double epsilon = 1e-5) {
    if (M1.size() != M2.size() || (M1.empty() ? 0 : M1[0].size()) != (M2.empty() ? 0 : M2[0].size())) {
        std::cerr << "Comparison failed: Matrix dimensions differ." << std::endl;
        return false;
    }
    if (M1.empty()) return true; // Both empty

    for (size_t i = 0; i < M1.size(); ++i) {
        for (size_t j = 0; j < M1[0].size(); ++j) {
            if (std::abs(M1[i][j] - M2[i][j]) > epsilon) {
                std::cerr << "Comparison failed at (" << i << "," << j << "): "
                          << M1[i][j] << " vs " << M2[i][j] << std::endl;
                return false;
            }
        }
    }
    return true;
}


// --- Formula Parsing and Application ---
struct ProductTerm {
    Matrix coeffs_A; // n x m
    Matrix coeffs_B; // m x p
    Matrix coeffs_C; // n x p (after transposing from file's p x n)
};

std::vector<double> parseRow(std::string row_str) {
    std::vector<double> row;
    // Remove '{' from start and '}' from end if they exist
    if (!row_str.empty() && row_str.front() == '{') row_str.erase(0, 1);
    if (!row_str.empty() && row_str.back() == '}') row_str.pop_back();

    std::stringstream ss(row_str);
    std::string item;
    while (std::getline(ss, item, ',')) {
        try {
            if (item.empty() && !ss.eof()) continue; // handle trailing comma before end or empty segments if format is loose
            if (item.empty() && ss.eof() && row.empty() && row_str.find_first_not_of(" \t\n\v\f\r") == std::string::npos) {
                 // This could be an empty row string like "{}" which becomes ""
                 // If it was meant to be an empty vector of doubles.
                 // For this specific format, rows are non-empty.
                 // throw std::runtime_error("Empty item in non-empty row string during parsing: '" + row_str + "'");
                 // If we expect all rows to have numbers, an empty item is an error or implies zero if that's the convention.
                 // The files usually have explicit zeros.
            }
            row.push_back(std::stod(item));
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error parsing number in row: '" << item << "' from '" << row_str << "'" << std::endl;
            throw;
        }  catch (const std::out_of_range& e) {
            std::cerr << "Error parsing number (out of range) in row: '" << item << "' from '" << row_str << "'" << std::endl;
            throw;
        }
    }
    return row;
}

Matrix parseCoeffMatrix(const std::string& mat_block_str) {
    Matrix mat;
    
    // Check for outer braces
    if (mat_block_str.length() < 4 || mat_block_str.substr(0, 2) != "{{" || mat_block_str.substr(mat_block_str.length() - 2) != "}}") {
        throw std::runtime_error("Invalid matrix block format (outer braces missing): " + mat_block_str);
    }
    
    // Parse the matrix by tracking brace levels
    size_t i = 1; // Start after the first '{'
    int matrix_brace_level = 1;
    
    while (i < mat_block_str.length() && matrix_brace_level > 0) {
        // Skip whitespace and commas
        while (i < mat_block_str.length() && (mat_block_str[i] == ' ' || mat_block_str[i] == ',')) {
            i++;
        }
        
        if (i >= mat_block_str.length()) break;
        
        // Check if we're at the end of the matrix
        if (mat_block_str[i] == '}') {
            matrix_brace_level--;
            if (matrix_brace_level == 0) break; // End of matrix
            i++;
            continue;
        }
        
        // We should be at the start of a row
        if (mat_block_str[i] != '{') {
            throw std::runtime_error("Expected '{' at start of row at position " + std::to_string(i) + ", found: " + std::string(1, mat_block_str[i]));
        }
        
        // Find the complete row
        size_t row_start = i;
        int row_brace_level = 0;
        while (i < mat_block_str.length()) {
            if (mat_block_str[i] == '{') {
                row_brace_level++;
            } else if (mat_block_str[i] == '}') {
                row_brace_level--;
                if (row_brace_level == 0) {
                    // Found the end of this row
                    std::string row_str = mat_block_str.substr(row_start, i - row_start + 1);
                    mat.push_back(parseRow(row_str));
                    i++; // Move past the closing brace
                    break;
                }
            }
            i++;
        }
        
        if (row_brace_level != 0) {
            throw std::runtime_error("Unmatched braces in row starting at position " + std::to_string(row_start));
        }
    }
    
    return mat;
}

std::vector<ProductTerm> loadFormula(const std::string& filepath, int n, int m, int p) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open formula file: " + filepath);
    }

    std::string raw_content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    file.close();

    // Basic cleanup: remove newlines, tabs, and spaces
    std::string content;
    content.reserve(raw_content.length());
    for(char c : raw_content) {
        if (c != '\n' && c != '\r' && c != '\t' && c != ' ') {
            content.push_back(c);
        }
    }

    std::vector<ProductTerm> formula;

    if (content.length() < 2 || content.front() != '{' || content.back() != '}') {
        throw std::runtime_error("Formula file content not enclosed in { }");
    }
    content = content.substr(1, content.length() - 2); // Remove top-level list braces
    if (content.empty()) return formula; // Empty list of terms

    size_t current_pos = 0;
    while (current_pos < content.length()) {
        if (content[current_pos] != '{') throw std::runtime_error("Expected '{' at start of term block");
        // Find the matching '}' for the current term's outer block {{{...}}}
        int brace_level = 0;
        size_t term_block_end = std::string::npos;
        for(size_t i = current_pos; i < content.length(); ++i) {
            if (content[i] == '{') brace_level++;
            else if (content[i] == '}') brace_level--;
            if (brace_level == 0) {
                term_block_end = i;
                break;
            }
        }
        if (term_block_end == std::string::npos) throw std::runtime_error("Unmatched braces for term block.");

        std::string term_str = content.substr(current_pos, term_block_end - current_pos + 1);
        current_pos = term_block_end + 1;
        if(current_pos < content.length() && content[current_pos] == ',') current_pos++; // Skip comma

        // term_str should be like {{{A_coeffs}}, {{B_coeffs}}, {{C_coeffs_file}}}
        // Note: The actual format has spaces after commas
        if (term_str.length() < 2 || term_str.front() != '{' || term_str.back() != '}') {
            throw std::runtime_error("Malformed term string: missing outer braces");
        }
        
        // Remove the outermost { and }
        std::string inner_str = term_str.substr(1, term_str.length() - 2);
        
        // Find the three matrices by counting braces properly
        int brace_count = 0;
        size_t matrix_start = 0;
        std::vector<std::string> matrices;
        
        for (size_t i = 0; i < inner_str.length(); ++i) {
            if (inner_str[i] == '{') {
                if (brace_count == 0) matrix_start = i;
                brace_count++;
            } else if (inner_str[i] == '}') {
                brace_count--;
                if (brace_count == 0) {
                    // Found the end of a matrix
                    matrices.push_back(inner_str.substr(matrix_start, i - matrix_start + 1));
                    if (matrices.size() == 3) break;
                }
            }
        }
        
        if (matrices.size() != 3) {
            throw std::runtime_error("Expected exactly 3 matrices in term, found " + std::to_string(matrices.size()));
        }
        
        std::string str_A = matrices[0];
        std::string str_B = matrices[1];
        std::string str_C_file = matrices[2];

        ProductTerm term;
        term.coeffs_A = parseCoeffMatrix(str_A);
        term.coeffs_B = parseCoeffMatrix(str_B);
        Matrix parsed_coeffs_C_from_file = parseCoeffMatrix(str_C_file);

        // Validate dimensions as they are parsed from file
        if (term.coeffs_A.size() != static_cast<size_t>(n) || (!term.coeffs_A.empty() && term.coeffs_A[0].size() != static_cast<size_t>(m)))
            throw std::runtime_error("coeffs_A dimensions mismatch (expected " + std::to_string(n) + "x" + std::to_string(m) + ")");
        if (term.coeffs_B.size() != static_cast<size_t>(m) || (!term.coeffs_B.empty() && term.coeffs_B[0].size() != static_cast<size_t>(p)))
            throw std::runtime_error("coeffs_B dimensions mismatch (expected " + std::to_string(m) + "x" + std::to_string(p) + ")");
        if (parsed_coeffs_C_from_file.size() != static_cast<size_t>(p) || (!parsed_coeffs_C_from_file.empty() && parsed_coeffs_C_from_file[0].size() != static_cast<size_t>(n)))
             throw std::runtime_error("parsed_coeffs_C_from_file dimensions mismatch (expected " + std::to_string(p) + "x" + std::to_string(n) +
                                     ", got " + std::to_string(parsed_coeffs_C_from_file.size()) + "x" + (parsed_coeffs_C_from_file.empty() ? "0" : std::to_string(parsed_coeffs_C_from_file[0].size())) + ")");


        // Transpose coeffs_C_from_file (p x n) to term.coeffs_C (n x p)
        term.coeffs_C = createMatrix(n, p);
        for (int r_idx = 0; r_idx < n; ++r_idx) {
            for (int c_idx = 0; c_idx < p; ++c_idx) {
                term.coeffs_C[r_idx][c_idx] = parsed_coeffs_C_from_file[c_idx][r_idx];
            }
        }
        formula.push_back(term);
    }

    std::cout << "Successfully parsed " << formula.size() << " product terms." << std::endl;
    return formula;
}

Matrix multiplyWithFormula(const Matrix& A, const Matrix& B, const std::vector<ProductTerm>& formula, int n_dim, int m_dim, int p_dim) {
    if (formula.empty()) {
         if (n_dim == 0 || p_dim == 0) return createMatrix(n_dim, p_dim); // Multiplying to/from zero-dim matrix
        throw std::runtime_error("Formula is empty for non-trivial multiplication.");
    }
    if (A.size() != static_cast<size_t>(n_dim) || (!A.empty() && A[0].size() != static_cast<size_t>(m_dim)))
        throw std::runtime_error("Input A dimensions mismatch provided n, m.");
    if (B.size() != static_cast<size_t>(m_dim) || (!B.empty() && B[0].size() != static_cast<size_t>(p_dim)))
        throw std::runtime_error("Input B dimensions mismatch provided m, p.");

    Matrix C = createMatrix(n_dim, p_dim);

    for (const auto& term : formula) {
        // Basic check: ensure formula terms are consistent with n_dim, m_dim, p_dim
        if (term.coeffs_A.size() != static_cast<size_t>(n_dim) || (!term.coeffs_A.empty() && term.coeffs_A[0].size() != static_cast<size_t>(m_dim)))
            throw std::runtime_error("Term coeffs_A dimensions mismatch during application.");
        if (term.coeffs_B.size() != static_cast<size_t>(m_dim) || (!term.coeffs_B.empty() && term.coeffs_B[0].size() != static_cast<size_t>(p_dim)))
            throw std::runtime_error("Term coeffs_B dimensions mismatch during application.");
        if (term.coeffs_C.size() != static_cast<size_t>(n_dim) || (!term.coeffs_C.empty() && term.coeffs_C[0].size() != static_cast<size_t>(p_dim)))
            throw std::runtime_error("Term coeffs_C dimensions mismatch during application (post-transpose).");

        double L_A = 0.0;
        for (int i = 0; i < n_dim; ++i) {
            for (int j = 0; j < m_dim; ++j) {
                if (term.coeffs_A[i][j] != 0) { // Optimization for sparse coeffs
                    L_A += term.coeffs_A[i][j] * A[i][j];
                }
            }
        }

        double L_B = 0.0;
        for (int i = 0; i < m_dim; ++i) {
            for (int j = 0; j < p_dim; ++j) {
                 if (term.coeffs_B[i][j] != 0) { // Optimization
                    L_B += term.coeffs_B[i][j] * B[i][j];
                }
            }
        }

        double P_k = L_A * L_B;

        if (P_k != 0) { // Optimization
            for (int i = 0; i < n_dim; ++i) {
                for (int j = 0; j < p_dim; ++j) {
                    if (term.coeffs_C[i][j] != 0) { // Optimization
                        C[i][j] += term.coeffs_C[i][j] * P_k;
                    }
                }
            }
        }
    }
    return C;
}


int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <formula_file.m> <n> <m> <p> [--verify]" << std::endl;
        std::cerr << "  n: rows of A, rows of C" << std::endl;
        std::cerr << "  m: cols of A, rows of B" << std::endl;
        std::cerr << "  p: cols of B, cols of C" << std::endl;
        std::cerr << "Example for 456: " << argv[0] << " matrix-multiplication/456/k0e7ba384e845ae2.m 4 5 6" << std::endl;
        return 1;
    }

    std::string formula_path = argv[1];
    int n = 0, m = 0, p = 0;
    bool verify = false;

    try {
        n = std::stoi(argv[2]);
        m = std::stoi(argv[3]);
        p = std::stoi(argv[4]);
        if (argc > 5 && std::string(argv[5]) == "--verify") {
            verify = true;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error parsing dimensions: " << e.what() << std::endl;
        return 1;
    }

    if (n <= 0 || m <= 0 || p <= 0) {
        std::cerr << "Dimensions must be positive." << std::endl;
        return 1;
    }

    std::cout << "Loading formula from: " << formula_path << std::endl;
    std::cout << "Multiplying " << n << "x" << m << " (A) by " << m << "x" << p << " (B) to get " << n << "x" << p << " (C)." << std::endl;

    std::vector<ProductTerm> formula;
    try {
        formula = loadFormula(formula_path, n, m, p);
    } catch (const std::exception& e) {
        std::cerr << "Error loading or parsing formula: " << e.what() << std::endl;
        return 1;
    }

    Matrix A = createMatrix(n, m);
    Matrix B = createMatrix(m, p);

    fillMatrixRandom(A);
    fillMatrixRandom(B);

    bool print_small_matrices = (n * m <= 30 && m * p <= 30 && n * p <= 30);

    if (print_small_matrices) {
        std::cout << "\nInput Matrix A:" << std::endl; printMatrix(A, "A");
        std::cout << "\nInput Matrix B:" << std::endl; printMatrix(B, "B");
    }


    Matrix C_formula;
    try {
        auto start_time = std::chrono::high_resolution_clock::now();
        C_formula = multiplyWithFormula(A, B, formula, n, m, p);
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = end_time - start_time;
        std::cout << "\nMultiplication with formula took: " << duration.count() << " ms (" << formula.size() << " terms)" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error during formula multiplication: " << e.what() << std::endl;
        return 1;
    }

    if (print_small_matrices) {
        std::cout << "\nResult Matrix C (from formula):" << std::endl; printMatrix(C_formula, "C_formula");
    }

    if (verify) {
        std::cout << "\n--- Verification with Naive Multiplication ---" << std::endl;
        try {
            auto start_naive_t = std::chrono::high_resolution_clock::now();
            Matrix C_naive = naiveMultiply(A, B);
            auto end_naive_t = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> naive_duration = end_naive_t - start_naive_t;
            std::cout << "Naive multiplication took: " << naive_duration.count() << " ms" << std::endl;

            if (print_small_matrices) {
                std::cout << "\nResult Matrix C (naive):" << std::endl; printMatrix(C_naive, "C_naive");
            }
            if (compareMatrices(C_formula, C_naive)) {
                std::cout << "Verification SUCCESSFUL: Formula result matches naive result." << std::endl;
            } else {
                std::cout << "Verification FAILED: Formula result DOES NOT match naive result." << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error during naive multiplication for verification: " << e.what() << std::endl;
        }
    }

    return 0;
}