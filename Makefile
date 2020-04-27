# Define a default compiler. This can be overridden via (say) 'make CXX=clang++'
CXX = g++-9
# Default flags are for users, not developers. Your code is not for you!
# -O3 turns on optimizations.
# -g places symbols in the binary.
#    It helps the users debug their code, and binary sizes aren't an issue these days.
# -march=native: Tells the compiler to generate the most sophisticated instructions for the CPU.
#    It is a good default for open-source software,
#    but maybe not for closed source, where binaries are delivered to clients.
# -Wfatal-errors: Stop compilation on first error.
#    Since C++ is well-nigh impossible to lex and parse even in the absence of syntax errors,
#    defaulting to compiling past errors is, in my opinion, insane.
#    But a default it is: If you forget this flag, prepare to be horrified by compiler messages.
# -ffast-math: Sadly, often necessary to get SIMD vectorization, but dangerous.
# -fno-finite-math-only: Disable bizarre optimizations that assume nans and infs never occur, which are enabled by -ffast-math
# -MMD: Emit header dependencies, but don't include system headers.
CXXFLAGS = -g -Wall -Wextra --std=gnu++17 -fno-finite-math-only -march=native -Wfatal-errors -MMD
LINKFLAGS = -L/usr/local/lib
INCFLAGS = -I/usr/local/include -I./lodepng -I./include

# Developers should compile with 'make DEBUG=1'
ifdef DEBUG
	# Compiling with sanitizers during development is simply essential:
	CXXFLAGS += -O2 -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer
else
	CXXFLAGS += -O3
endif

# Find all .cpp files:
SRCS := $(wildcard *.cpp)
# For each .cpp, create a foo.x for foo.cpp:
EXECS := $(patsubst %.cpp,%.x,$(SRCS))

all: $(EXECS)

#This builds all .cpp files into separate executables:
%.x: %.cpp
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -o $@ $< lodepng/lodepng.cpp $(LINKFLAGS)

-include $(SRCS:.cpp=.d)

.PHONY: clean
clean:
	rm -rf *.png *.x *.x.dSYM/ *.d examples/*.x
