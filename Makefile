CXX      := g++
INC      := -I ./include

# Test suite settings
BUILD       := build
HEADERS     := $(shell find include -name '*.hpp')
TEST_SRC    := $(wildcard tests/*.cpp)
TEST_DEPS   := $(TEST_SRC) tests/matl_test.hpp $(HEADERS)
TEST_FLAGS  := -std=c++20 -O1 -g -Wall -Wextra -Werror
SAN_FLAGS   := -fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer

# Benchmark settings (optimized, without MATL debug checks and trace output)
BENCH_SRC   := $(wildcard benchmarks/*.cpp)
BENCH_FLAGS := -std=c++20 -O2 -Wall -Wextra -Werror -DNDEBUG

# Example programs (examples/*.cpp, built into build/examples/<name>, with MATL's debug checks on)
EXAMPLE_SRC   := $(wildcard examples/*.cpp)
EXAMPLE_BIN   := $(patsubst examples/%.cpp,$(BUILD)/examples/%,$(EXAMPLE_SRC))
EXAMPLE_FLAGS := -std=c++20 -O1 -g -Wall -Wextra -Werror

.PHONY: test clean examples run-examples tests check check-sanitize check-nodebug check-abort check-noexceptions check-openmp benchmarks bench

# Run the test suite (same as 'make check')
test: check

# Build the example programs
examples: $(EXAMPLE_BIN)

$(BUILD)/examples/%: examples/%.cpp $(HEADERS)
	@mkdir -p $(dir $@)
	$(CXX) $(EXAMPLE_FLAGS) $< $(INC) -o $@

# Build and run the example programs (their output is discarded; a failing example stops the run)
run-examples: examples
	@for e in $(EXAMPLE_BIN); do ./$$e > /dev/null && echo "ok      $$e" || { echo "FAILED  $$e"; exit 1; }; done

# Build the test suite
tests: $(BUILD)/matl_tests

$(BUILD)/matl_tests: $(TEST_DEPS)
	@mkdir -p $(BUILD)
	$(CXX) $(TEST_FLAGS) $(TEST_SRC) $(INC) -I tests -o $@

$(BUILD)/matl_tests_sanitize: $(TEST_DEPS)
	@mkdir -p $(BUILD)
	$(CXX) $(TEST_FLAGS) $(SAN_FLAGS) -DMATL_DEBUG $(TEST_SRC) $(INC) -I tests -o $@

# Run the test suite. Results go to stderr; anything the library prints to stdout goes to a log file.
# Pass a filter to run a subset of test cases, e.g. 'make check FILTER=Transpose'
check: $(BUILD)/matl_tests
	./$(BUILD)/matl_tests $(FILTER) > $(BUILD)/matl_tests.stdout.log

# Run the test suite with MATL_NO_DEBUG (no bounds/dimension checks or trace output)
$(BUILD)/matl_tests_nodebug: $(TEST_DEPS)
	@mkdir -p $(BUILD)
	$(CXX) $(TEST_FLAGS) -DMATL_NO_DEBUG $(TEST_SRC) $(INC) -I tests -o $@

check-nodebug: $(BUILD)/matl_tests_nodebug
	./$(BUILD)/matl_tests_nodebug $(FILTER)

# Run the test suite with MATL_ABORT_ON_ERROR (errors abort instead of throwing)
$(BUILD)/matl_tests_abort: $(TEST_DEPS)
	@mkdir -p $(BUILD)
	$(CXX) $(TEST_FLAGS) -DMATL_ABORT_ON_ERROR $(TEST_SRC) $(INC) -I tests -o $@

check-abort: $(BUILD)/matl_tests_abort
	./$(BUILD)/matl_tests_abort $(FILTER)

# Run the test suite with OpenMP (-fopenmp; GCC, or Clang with an OpenMP runtime), on 4 threads
$(BUILD)/matl_tests_openmp: $(TEST_DEPS)
	@mkdir -p $(BUILD)
	$(CXX) $(TEST_FLAGS) -fopenmp $(TEST_SRC) $(INC) -I tests -o $@

check-openmp: $(BUILD)/matl_tests_openmp
	OMP_NUM_THREADS=4 ./$(BUILD)/matl_tests_openmp $(FILTER) > $(BUILD)/matl_tests_openmp.stdout.log

# Check that MATL compiles without exceptions, and that errors then abort
$(BUILD)/no_exceptions: tests/no_exceptions/main.cpp $(HEADERS)
	@mkdir -p $(BUILD)
	$(CXX) $(TEST_FLAGS) -fno-exceptions $< $(INC) -o $@

check-noexceptions: $(BUILD)/no_exceptions
	./$(BUILD)/no_exceptions

# Run the test suite with AddressSanitizer and UndefinedBehaviorSanitizer (and the MATL_DEBUG
# trace enabled, so that the trace code is compiled and exercised; its output goes to a log file).
# alloc_dealloc_mismatch is on by default on Linux but not on macOS, so it is set here.
check-sanitize: $(BUILD)/matl_tests_sanitize
	ASAN_OPTIONS=alloc_dealloc_mismatch=1 ./$(BUILD)/matl_tests_sanitize $(FILTER) > $(BUILD)/matl_tests_sanitize.stdout.log

# Build the benchmarks
benchmarks: $(patsubst benchmarks/%.cpp,$(BUILD)/benchmarks/%,$(BENCH_SRC))

$(BUILD)/benchmarks/%: benchmarks/%.cpp $(HEADERS)
	@mkdir -p $(dir $@)
	$(CXX) $(BENCH_FLAGS) $< $(INC) -o $@

# Build and run the benchmarks
bench: benchmarks
	@for b in $(patsubst benchmarks/%.cpp,$(BUILD)/benchmarks/%,$(BENCH_SRC)); do ./$$b; done

clean:
	rm -rf $(BUILD)
