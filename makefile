all:	rm_reads.cpp
	c++ -std=c++0x -g -O0 -D DEBUG rm_reads.cpp -o rm_reads
clean:	
	rm rm_reads
