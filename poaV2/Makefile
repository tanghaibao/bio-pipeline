
AR=ar rc

TARGETS=poa liblpo.a poa_doc libbflag.a

# align_score.c CAN BE USED TO ADD CUSTOMIZED SCORING FUNCTIONS
OBJECTS= \
	align_score.o \
	main.o


LIBOBJECTS= \
	black_flag.o \
	seq_util.o \
	fasta_format.o \
	msa_format.o \
	align_lpo2.o \
	align_lpo_po2.o \
	buildup_lpo.o \
	lpo.o \
	heaviest_bundle.o \
	lpo_format.o \
	create_seq.o \
	remove_bundle.o \
	numeric_data.o \
	stringptr.o


CC = gcc
#CFLAGS= -g -ansi-strict -W -Wall -DUSE_WEIGHTED_LINKS -DUSE_PROJECT_HEADER -I.
CFLAGS= -g -DUSE_WEIGHTED_LINKS -DUSE_PROJECT_HEADER -I.
# -I$(HOME)/lib/include
# -DREPORT_MAX_ALLOC

clean:
	rm -f $(OBJECTS) $(LIBOBJECTS) $(TARGETS)

liblpo.a: $(LIBOBJECTS)
	rm -f $@
	$(AR) $@ $(LIBOBJECTS)
	ranlib $@



# NB: LIBRARY MUST FOLLOW OBJECTS OR LINK FAILS WITH UNRESOLVED REFERENCES!!
poa: $(OBJECTS) liblpo.a
	$(CC) -o $@ $(OBJECTS) -lm liblpo.a

what:
	@echo poa: partial-order based sequence alignment program
	@echo liblpo.a: partial-order alignment and utilities function library
