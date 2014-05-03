CXX= g++
CXXFLAGS = -std=c++0x -Wall
OPT = -O2
DEBUG = -g -O0 -D DEBUG
all: src/rm_reads.cpp
	$(CXX) $(CXXFLAGS) $(OPT) src/*.cpp -o rm_reads
.PHONY: clean
clean:
	rm -rf rm_reads

