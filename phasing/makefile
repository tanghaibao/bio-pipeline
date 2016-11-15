OS := $(shell uname)
ifeq ($(OS), Darwin)
	LDFLAGS+=-lpthread -lz -lcurl -lcrypto
else
	LDFLAGS+=-lpthread -lrt -lz -lcurl -lcrypto
endif

CXX=g++-6
CXXFLAGS+=-Ihtslib -Ihtslib/htslib -Iseqan/include -std=c++14 -DNDEBUG -W -Wall -pedantic -fopenmp

default: all
all: phasing

htslib/libhts.a:
	cd htslib && autoheader && autoconf && ./configure --enable-libcurl && $(MAKE) -j 32

phasing: main.o htslib/libhts.a
	$(CXX) -o $@ $^ $(LDFLAGS)

main.o: main.cc
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

clean:
	rm -f main.o main

.PHONY: default all clean
