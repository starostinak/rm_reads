CXX= g++
CXXFLAGS = -std=c++0x -Wall
DEBUG = -g -O0 -D DEBUG
all: *.cpp
	$(CXX) $(CXXFLAGS) *.cpp -o rm_reads
.PHONY: clean
clean:
	rm -rf rm_reads

