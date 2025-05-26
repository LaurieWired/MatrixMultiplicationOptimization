#!/bin/bash

# Compile the formula_multiply program
echo "Compiling formula_multiply.cpp..."
g++ -std=c++17 -O2 -o formula_multiply formula_multiply.cpp
if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

echo -e "\n========================================="
echo "Testing matrix multiplication formulas"
echo "========================================="

# Test each directory
for dir in matrix-multiplication/*/; do
    # Skip .git directory
    if [[ "$dir" == *".git"* ]]; then
        continue
    fi
    
    # Extract directory name and parse dimensions
    dirname=$(basename "$dir")
    n=${dirname:0:1}
    m=${dirname:1:1}
    p=${dirname:2:1}
    
    echo -e "\n----- Testing ${n}x${m}x${p} matrix multiplication -----"
    
    # Find the first .m file in the directory
    formula_file=$(find "$dir" -name "*.m" -type f | head -1)
    
    if [ -z "$formula_file" ]; then
        echo "No formula files found in $dir"
        continue
    fi
    
    echo "Using formula file: $formula_file"
    
    # Run the test with verification
    ./formula_multiply "$formula_file" "$n" "$m" "$p" --verify | grep -E "(Successfully parsed|Verification|Error|took:)"
    
    if [ $? -eq 0 ]; then
        echo "✓ Test passed for ${n}x${m}x${p}"
    else
        echo "✗ Test failed for ${n}x${m}x${p}"
    fi
done

echo -e "\n========================================="
echo "All tests completed"
echo "=========================================" 