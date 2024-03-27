# Define compiler and flags
CXX = g++
CXXFLAGS = -std=c++20 -ggdb -fopenmp -pedantic -Wall -Weffc++ -Wextra -Wconversion -Wsign-conversion

# Define source and object files
SRC = main.cpp
OBJ = $(SRC:.cpp=.o)

# Default target - builds the final executable
output: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o output

# Compile source code into object files
$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Remove object files for cleaning
clean:
	rm -f $(OBJ)
