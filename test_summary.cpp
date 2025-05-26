#include <iostream>
#include <iomanip>
#include <vector>

struct TestResult {
    int n, m, p;
    int formula_terms;
    int naive_ops;
    double efficiency_gain;
};

int main() {
    std::vector<TestResult> results = {
        {4, 5, 6, 90, 4*5*6, (4*5*6 - 90) / double(4*5*6) * 100},
        {4, 5, 7, 104, 4*5*7, (4*5*7 - 104) / double(4*5*7) * 100},
        {4, 6, 6, 106, 4*6*6, (4*6*6 - 106) / double(4*6*6) * 100},
        {4, 6, 7, 123, 4*6*7, (4*6*7 - 123) / double(4*6*7) * 100},
        {4, 7, 7, 144, 4*7*7, (4*7*7 - 144) / double(4*7*7) * 100},
        {5, 5, 6, 110, 5*5*6, (5*5*6 - 110) / double(5*5*6) * 100},
        {5, 5, 7, 127, 5*5*7, (5*5*7 - 127) / double(5*5*7) * 100},
        {5, 6, 6, 130, 5*6*6, (5*6*6 - 130) / double(5*6*6) * 100},
        {5, 6, 7, 150, 5*6*7, (5*6*7 - 150) / double(5*6*7) * 100},
        {5, 7, 7, 176, 5*7*7, (5*7*7 - 176) / double(5*7*7) * 100},
        {6, 6, 7, 183, 6*6*7, (6*6*7 - 183) / double(6*6*7) * 100}
    };
    
    std::cout << "\n============================================" << std::endl;
    std::cout << "Matrix Multiplication Efficiency Summary" << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << std::setw(10) << "Size" 
              << std::setw(15) << "Naive Ops"
              << std::setw(15) << "Formula Ops"
              << std::setw(15) << "Reduction %" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    
    for (const auto& r : results) {
        std::cout << std::setw(8) << (std::to_string(r.n) + "x" + std::to_string(r.m) + "x" + std::to_string(r.p))
                  << std::setw(15) << r.naive_ops
                  << std::setw(15) << r.formula_terms
                  << std::setw(14) << std::fixed << std::setprecision(1) << r.efficiency_gain << "%" << std::endl;
    }
    
    std::cout << "\nAll formulas successfully reduce the number of" << std::endl;
    std::cout << "scalar multiplications compared to naive algorithm!" << std::endl;
    
    return 0;
} 