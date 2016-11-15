OS := $(shell uname)
ifeq ($(OS), Darwin)
	LDFLAGS+=-lpthread -lz -lcurl -lcrypto
else
	LDFLAGS+=-lpthread -lrt -lz -lcurl -lcrypto
endif

PROG=phasing
CXX=g++-6
CXXFLAGS+=-Ihtslib -Ihtslib/htslib -Iseqan/include -std=c++14 -DNDEBUG -W -Wall -pedantic -fopenmp
SRCS=main.cc program.cc
OBJS=$(SRCS:.cc=.o)

default: all
all: phasing

htslib/libhts.a:
	cd htslib && autoheader && autoconf && ./configure --enable-libcurl && $(MAKE) -j 32

$(PROG): $(OBJS) htslib/libhts.a
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cc
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

clean:
	@rm -f $(PROG) *.o

.PHONY: default all clean
