# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -Wno-sign-compare -Wno-endif-labels -Wno-unused-variable -O3 -std=c++14

# Include directories for GSL and other headers
INCLUDES = -I/usr/include

# Source files
SRCS = constant.cpp dynamics.cpp debye_bath.cpp ohmic_bath.cpp polar.cpp TD_hamiltonian.cpp param.cpp super_ohmic_bath.cpp utilize.cpp

# Header files
HDRS = constant.h dynamics.h debye_bath.h ohmic_bath.h polar.h TD_hamiltonian.h param.h super_ohmic_bath.h utilize.h

# Object files
OBJS = $(SRCS:.cpp=.o)

# Static library name
LIB = libcpu2des.a

# Default target
all: $(LIB)

# Compile object files
%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Create static library
$(LIB): $(OBJS)
	ar rcs $(LIB) $(OBJS)

# Clean up the build directory
clean:
	rm -f $(OBJS) $(LIB) *.a

# Phony targets
.PHONY: all clean