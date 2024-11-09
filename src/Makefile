# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -Wno-sign-compare -Wno-endif-labels -Wno-unused-variable -O3 -std=c++14

# Include directories for GSL, other headers, and the library directory. change as need
INCLUDES = -I/usr/include -I../lib

# Libraries to link against
LIBS = -L../lib -lgsl -lgslcblas -lm -lcpu2des

# Source files
SRCS = CPU_2DES.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
EXEC = CPU_2DES

# Default target
all: $(EXEC)

# Link the executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJS) $(LIBS)

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean up the build directory
clean:
	rm -f $(OBJS) $(EXEC)

# Phony targets
.PHONY: all clean
