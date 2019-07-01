CC=gcc
# PROF=-pg
CFLAGS=-Wall -g $(PROF) -O2 -I.
LDFLAGS=-g $(PROF)

.PHONY : all test clean 

all: sqr_test sqr

sqrlib.o : sqrlib.c
	$(CC) -c $(CFLAGS) -o $@ $<

sqr.o : sqr.c
	$(CC) -c $(CFLAGS) -o $@ $<

sqr_test.o : test/sqr_test.c
	$(CC) -c $(CFLAGS) -o $@ $<

cpuid.o : test/cpuid.c
	$(CC) -c $(CFLAGS) -o $@ $<

sqr_test: sqrlib.o sqr_test.o cpuid.o
	$(CC) -lcrypto -o $@ $^ ${LDFLAGS}

sqr: sqrlib.o sqr.o
	$(CC) -lcrypto -o $@ $^ ${LDFLAGS}

test: sqr_test
	./sqr_test

valgrind: sqr_test
	valgrind --leak-check=full ./sqr_test

style::
	uncrustify -c style/uncrustify.cfg -f sqrlib.c -o sqrlib.c
	uncrustify -c style/uncrustify.cfg -f test/sqr_test.c -o test/sqr_test.c
	uncrustify -c style/uncrustify.cfg -f test/cpuid.c -o test/cpuid.c

clean:
	-rm -f *.o sqr_test sqr
