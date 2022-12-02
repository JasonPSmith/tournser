build: tournser

tournser: tournser.cpp
	@echo "Compiling \"tournser\"." && g++ -std=c++14 -O3 -pthread tournser.cpp -o tournser

clean:
	rm -f tournser
