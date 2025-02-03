CXX = g++
CXXFLAGS = -std=c++20 -I/usr/include/eigen3 -I/usr/include -lboost_system -lboost_filesystem -lboost_math_c99 -Wall -Wextra
TARGET = output
SOURCES = main.cpp

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

clean:
	rm -f $(TARGET)
