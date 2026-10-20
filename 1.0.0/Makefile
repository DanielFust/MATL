CXX      := g++
CPPFLAGS := -std=c++20 -O0 -Wall -Wextra -g3
#CPPFLAGS := -std=c++20 -O0 -w -g3
INC      := -I ./include

.PHONY: test clean

test: 
	$(CXX) $(CPPFLAGS) ./examples/test.cpp $(INC) -o $@

clean:
	rm test
