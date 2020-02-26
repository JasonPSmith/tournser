build: tournser

tournser: tournser.cpp
	c++ -pthread -std=c++14 tournser.cpp -o tournser -Ofast -D NDEBUG
