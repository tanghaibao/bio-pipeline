PROGS=gbs_digest
SRCS=gbs_digest.c cmdline.c
OBJS=$(SRCS:.c=.o)
CC=gcc
LDFLAGS=-static -lz
CFLAGS=-O3 -Wall -I.

all: $(PROGS)

cmdline.h cmdline.c: gbs_digest.ggo
	gengetopt -i $< -a args_info -f parser -u

${PROGS}: $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

%.o: %.c cmdline.c
	$(CC) $(CFLAGS) $< -c -o $@

clean:
	rm -f cmdline.?
	rm -f $(PROGS)
	rm -f *.o

test: $(PROGS)
	./${PROGS} test.fq.gz
