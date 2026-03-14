# Variables - Define your "tuning knobs" here
CXX      = clang++
CXXFLAGS = -O3 -std=c++17 -Iinclude -Wall -Wextra
TARGET   = out/main

# Automatically find all .cpp files in the src directory
SRCS     = $(wildcard src/*.cpp)

# Translate 'src/filename.cpp' into 'out/filename.o'
OBJS     = $(patsubst src/%.cpp, out/%.o, $(SRCS))

# Default target
all: $(TARGET)

# The default rule: Build the final executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Rule for how to compile a single .cpp file into a .o file inside out/
# The '| out' is an order-only prerequisite, ensuring the folder exists first.
out/%.o: src/%.cpp | out
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to create the out folder if it doesn't exist
out:
	mkdir -p out

# Clean rule to delete the compiled files and start fresh
clean:
	rm -rf out/*.o $(TARGET)

# Phony targets (not actual files)
.PHONY: all clean
