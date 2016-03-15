CC = gcc
CFLAGS = -g -std=c99 -w -DDEBUG

TARGETS = tests_rbt tests tests_sparse_compose tests_dag tests_striping tests_mbr tests_msr

INCLUDE = ../*.c -I.. -I../../jerasure-kmg/src
LIB = ../../jerasure-kmg/src/*.o

all: $(TARGETS)

tests_rbt: tests_rbt.c test_common.h ../*.c
	$(CC) $(CFLAGS) $(LIB) tests_rbt.c $(INCLUDE) -o tests_rbt

tests_dag: tests_dag.c ../*.c
	$(CC) $(CFLAGS) $(LIB) tests_dag.c $(INCLUDE) -o tests_dag

tests: tests.c
	$(CC) $(CFLAGS) $(LIB) tests.c $(INCLUDE) -o tests

tests_sparse_compose: tests_sparse_compose.c
	$(CC) $(CFLAGS) $(LIB) tests_sparse_compose.c $(INCLUDE) -o tests_sparse_compose

tests_striping: tests_striping.c
	$(CC) $(CFLAGS) $(LIB) tests_striping.c $(INCLUDE) -o tests_striping

tests_mbr: tests_mbr.c
	$(CC) $(CFLAGS) $(LIB) tests_mbr.c $(INCLUDE) -o tests_mbr

tests_msr: tests_msr.c
	$(CC) $(CFLAGS) $(LIB) tests_msr.c $(INCLUDE) -o tests_msr

clean:
	rm -f $(TARGETS) *.o
	rm -rf *.dSYM
