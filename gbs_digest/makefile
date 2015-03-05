PROGS=gbs_digest
SRCS=gbs_digest.c cmdline.c
OBJS=$(SRCS:.c=.o)
CC=gcc
LDFLAGS=-static -lz
CFLAGS=-O3 -Wall -I.

all: $(PROGS)

cmdline.c: gbs_digest.ggo
	gengetopt -i $< -a args_info -f parser --unamed-opts

${PROGS}: $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) $< -c -o $@

clean:
	rm -f $(PROGS)
	rm -f *.o

test: $(PROGS)
	./${PROGS} "GATC|[AG](CATG)[CT]" test.fq.gz
