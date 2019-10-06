CC=gcc
# PROF=-pg
CFLAGS=-Wall -g $(PROF) -O2 -I.
CFLAGS_OPT=-Wall -g $(PROF) -O3 -I.
LDFLAGS=-g $(PROF)

V8_NAME=sqrlib_v8_avx512

.PHONY : all test clean 

all: sqr_test sqr

sqrlib.o : sqrlib/sqrlib.c
	$(CC) -c $(CFLAGS) -o $@ $<

$(V8_NAME).o : sqrlib/$(V8_NAME).c
	$(CC) -c $(CFLAGS_OPT) -mavx512ifma -mavx512bw -o $@ $<
	objdump -l -d -M att -S $@ > $@.asm

sqrlib.a: sqrlib.o $(V8_NAME).o
	ar -crv $@ $^

sqr.o : sqr.c
	$(CC) -c $(CFLAGS) -o $@ $<

sqr_test.o : test/sqr_test.c
	$(CC) -c $(CFLAGS) -o $@ $<

cpuid.o : test/cpuid.c
	$(CC) -c $(CFLAGS) -o $@ $<

sqr_test: sqr_test.o cpuid.o sqrlib.a
	$(CC) -lcrypto -o $@ $^ ${LDFLAGS}

sqr: sqr.o sqrlib.a
	$(CC) -lcrypto -o $@ $^ ${LDFLAGS}

test: sqr_test
	./sqr_test

test_sda: sqr_test
	~/bin/sde-external-8.35.0-2019-03-11-lin/sde64 -- ./sqr_test

valgrind: sqr_test
	valgrind --leak-check=full ./sqr_test

STYLE_OPT=--no-backup -c style/uncrustify.cfg

style::
	uncrustify $(STYLE_OPT) -f sqrlib/sqrlib.c -o sqrlib/sqrlib.c
	uncrustify $(STYLE_OPT) -f sqrlib/sqrlib_v8_avx512.c -o sqrlib/sqrlib_v8_avx512.c
	uncrustify $(STYLE_OPT) -f sqrlib/sqrlib_v8_avx512_modN.h -o sqrlib/sqrlib_v8_avx512_modN.h
	uncrustify $(STYLE_OPT) -f test/sqr_test.c -o test/sqr_test.c
	uncrustify $(STYLE_OPT) -f test/cpuid.c -o test/cpuid.c

clean:
	-rm -f *.o *.o.asm *.a sqr_test sqr
