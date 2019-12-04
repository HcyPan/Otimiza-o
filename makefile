CXX=clang++
CXXFLAGS=-std=c++11 -stdlib=libc++ -I headers -Wall -pedantic
BIN=bin

SRC=$(wildcard src/*.cpp)
OBJ=$(SRC:%.cpp=%.o)

all: build run 
		
build: $(OBJ)
		$(CXX) $(CXXFLAGS) -o $(BIN) main.cpp $^
		

%.o: %.c
		$(CXX) $@ -c $<

run: 
		./bin

clean:
		rm -f *.o