PROGS=gbs_digest
SRCS=gbs_digest.c
OBJS=$(SRCS:.c=.o)
CC=gcc
LDFLAGS=-static -lz
CFLAGS=-O3 -Wall -I.

all: $(PROGS)

${PROGS}: $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(PROGS)
	rm -f *.o

test: $(PROGS)
	./${PROGS} "GATC|[AG](CATG)[CT]" test.fq.gz
