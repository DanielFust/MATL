CXX      := g++
CPPFLAGS := -std=c++20 -O0 -Wall -Wextra -g3
#CPPFLAGS := -std=c++20 -O0 -w -g3
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

.PHONY: test clean examples tests check check-sanitize check-nodebug benchmarks bench

# Build the example program
test:
	$(CXX) $(CPPFLAGS) ./examples/test.cpp $(INC) -o $@

# Build the example program into the build directory
examples: $(BUILD)/examples/test

$(BUILD)/examples/test: examples/test.cpp $(HEADERS)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $< $(INC) -o $@

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

# Run the test suite with AddressSanitizer and UndefinedBehaviorSanitizer (and the MATL_DEBUG
# trace enabled, so that the trace code is compiled and exercised; its output goes to a log file)
check-sanitize: $(BUILD)/matl_tests_sanitize
	./$(BUILD)/matl_tests_sanitize $(FILTER) > $(BUILD)/matl_tests_sanitize.stdout.log

# Build the benchmarks
benchmarks: $(patsubst benchmarks/%.cpp,$(BUILD)/benchmarks/%,$(BENCH_SRC))

$(BUILD)/benchmarks/%: benchmarks/%.cpp $(HEADERS)
	@mkdir -p $(dir $@)
	$(CXX) $(BENCH_FLAGS) $< $(INC) -o $@

# Build and run the benchmarks
bench: benchmarks
	@for b in $(patsubst benchmarks/%.cpp,$(BUILD)/benchmarks/%,$(BENCH_SRC)); do ./$$b; done

clean:
	rm -f test
	rm -rf $(BUILD)
