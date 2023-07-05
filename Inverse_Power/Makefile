output: main.o
	g++ main.o -o output
    
main.o: main.cpp
	g++ -std=c++20 -ggdb -fopenmp -pedantic -c -Wall -Weffc++ -Wextra -Wconversion -Wsign-conversion main.cpp
    
clean:
	rm *.o
