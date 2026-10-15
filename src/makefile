CC      ?= cc
CFLAGS  ?= -O2 -Wall -pedantic
RM 	= rm
PREFIX  ?= /usr/local
BINDIR  = $(PREFIX)/bin
# SRC	= src/
# PY	= python

default: all

# .PHONY: default all concatenate subsampling collapser translate install run clean
.PHONY: default all concatenate collapser translate install run clean

# all: concatenate subsampling collapser translate
all: concatenate collapser translate

concatenate: concatenate_tables.c ht_array_DYNAMIC.c
	$(CC) $(CFLAGS) -o concatenate_tables.o concatenate_tables.c ht_array_DYNAMIC.c -lm

# subsampling: subsampling_seq.c getlines.c ht.c
# 	$(CC) $(CFLAGS) -o subsampling_seq.o subsampling_seq.c getlines.c ht.c -lm

collapser: hash_collapser.c list_counts.c ht.c
	$(CC) $(CFLAGS) -o hash_collapser.o hash_collapser.c list_counts.c ht.c -lm

translate: translate.c
	$(CC) $(CFLAGS) -o translate.o translate.c -lz

# Installs the executables and the bash scripts under $(PREFIX)/bin
install: all
	install -d $(BINDIR)
	install -m 755 concatenate_tables.o    $(BINDIR)/concatenate_tables
	# install -m 755 subsampling_seq.o $(BINDIR)/subsampling_seq
	install -m 755 hash_collapser.o   $(BINDIR)/hash_collapser
	install -m 755 translate.o $(BINDIR)/translate
	install -m 755 ../fast_trimming.sh ../comp_dna.sh $(BINDIR)




# Files and output name used by `make run`, override with e.g.
# $ make run RUN_FILES="a.counts b.counts c.counts" RUN_OUTPUT=my_output
RUN_FILES  ?= data_test/T4.counts data_test/T6.counts
RUN_OUTPUT ?= output_test

run: concatenate
	./concatenate_tables.o -o $(RUN_OUTPUT) $(RUN_FILES)

# test:
# 	$(PY) validate.py $(RUN_OUTPUT).csv $(RUN_FILES)

clean:
	$(RM) *.o
