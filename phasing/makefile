CXX=g++-6
CXXFLAGS+=-Iseqan/include -std=c++14 -DNDEBUG -O3 -W -Wall -pedantic -fopenmp
LDFLAGS+=-lpthread -lrt

default: all
all: main

main: main.o
	$(CXX) $(LDFLAGS) -o $@ $^

main.o: main.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $^

clean:
	rm -f main.o main

.PHONY: default all clean
