all:
	  g++ -g -pedantic -Wextra -Wshadow -Wnon-virtual-dtor -Wall -pthread *.cpp -o raytracer -std=c++11 -O3 

