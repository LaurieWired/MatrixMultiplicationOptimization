# Makefile for Matrix Multiplication Optimization Project

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra
LDFLAGS = 

# Source files
SOURCES = formula_multiply.cpp strassen_multiply.cpp test_summary.cpp naive_multiply.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLES = formula_multiply strassen_multiply test_summary naive_multiply

# Default target
all: $(EXECUTABLES)

# Build formula_multiply
formula_multiply: formula_multiply.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Build strassen_multiply
strassen_multiply: strassen_multiply.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Build test_summary
test_summary: test_summary.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Build naive_multiply
naive_multiply: naive_multiply.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Generic rule for object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Run tests
test: formula_multiply
	@echo "Running tests on all matrix sizes..."
	@./test_all_sizes.sh

# Run summary
summary: test_summary
	@./test_summary

# Run Strassen benchmark
benchmark: strassen_multiply
	@echo "Running Strassen vs Naive benchmark..."
	@./strassen_multiply

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(EXECUTABLES)

# Install (optional - you can customize the install path)
install: $(EXECUTABLES)
	@echo "Installing executables to /usr/local/bin (requires sudo)"
	@echo "Run 'sudo make install' if you want to install system-wide"
	install -m 755 $(EXECUTABLES) /usr/local/bin/

# Uninstall
uninstall:
	rm -f /usr/local/bin/formula_multiply /usr/local/bin/strassen_multiply /usr/local/bin/test_summary

# Help target
help:
	@echo "Available targets:"
	@echo "  all          - Build all executables (default)"
	@echo "  formula_multiply - Build the formula multiply program"
	@echo "  strassen_multiply - Build the Strassen algorithm benchmark"
	@echo "  naive_multiply - Build the naive multiplication program"
	@echo "  test_summary - Build the test summary program"
	@echo "  test         - Run tests on all matrix sizes"
	@echo "  summary      - Display efficiency summary"
	@echo "  benchmark    - Run Strassen vs Naive benchmark"
	@echo "  clean        - Remove build artifacts"
	@echo "  install      - Install executables to /usr/local/bin"
	@echo "  uninstall    - Remove installed executables"
	@echo "  help         - Display this help message"

# Mark targets that don't create files
.PHONY: all test summary benchmark clean install uninstall help 