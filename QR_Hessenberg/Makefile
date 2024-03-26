# Makefile for Writing Make Files Example
 
# *****************************************************
# Variables to control Makefile operation
 
CC = g++
CFLAGS = -std=c++20 -ggdb -fopenmp -pedantic -Wall -Weffc++ -Wextra -Wconversion -Wsign-conversion
 
# ****************************************************ru
# Targets needed to bring the executable up to date
 
output: main.o
	$(CC) $(CFLAGS) -o output main.o
 
# The main.o target can be written more simply
 
main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

clean:
	rm *.o

run: 
	./output
