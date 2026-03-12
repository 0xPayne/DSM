# Variables - Define your "tuning knobs" here
CXX      = clang++
CXXFLAGS = -O3 -std=c++17 -Iinclude -Wall -Wextra
TARGET   = out/main

# Automatically find all .cpp files in the src directory
SRCS     = $(wildcard src/*.cpp)

# Define the object files (.o) based on the .cpp files
OBJS     = $(SRCS:.cpp=.o)

# The default rule: Build the final executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Rule for how to compile a single .cpp file into a .o file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule to delete the compiled files and start fresh
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets (not actual files)
.PHONY: all clean
