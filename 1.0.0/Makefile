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

.PHONY: test clean examples tests check check-sanitize

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
	$(CXX) $(TEST_FLAGS) $(SAN_FLAGS) $(TEST_SRC) $(INC) -I tests -o $@

# Run the test suite. Results go to stderr; MATL_DEBUG trace output (stdout) goes to a log file.
# Pass a filter to run a subset of test cases, e.g. 'make check FILTER=Transpose'
check: $(BUILD)/matl_tests
	./$(BUILD)/matl_tests $(FILTER) > $(BUILD)/matl_tests.stdout.log

# Run the test suite with AddressSanitizer and UndefinedBehaviorSanitizer
check-sanitize: $(BUILD)/matl_tests_sanitize
	./$(BUILD)/matl_tests_sanitize $(FILTER) > $(BUILD)/matl_tests_sanitize.stdout.log

clean:
	rm -f test
	rm -rf $(BUILD)
