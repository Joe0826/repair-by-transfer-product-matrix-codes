CC = gcc
CFLAGS = -O3 -std=c99

TARGETS = encoder_pmc collector_pmc repair_pmc tests encoder_blk_pmc

LIB = ../*.o *.o -lgf_complete -lJerasure -lm -Wl,-rpath -Wl,/usr/local/lib
#../../gf-complete-kmg/*.o
#../../jerasure-kmg/*.o
INCLUDE = -I.. -I../../jerasure-kmg/include


all: clean $(TARGETS)

encoder_pmc: encoder_pmc.c common
	$(CC) $(CFLAGS) $(LIB) encoder_pmc.c $(INCLUDE) -o encoder_pmc

encoder_blk_pmc: encoder_blk_pmc.c common
	$(CC) $(CFLAGS) $(LIB) encoder_blk_pmc.c $(INCLUDE) -o encoder_blk_pmc

collector_pmc: collector_pmc.c common
	$(CC) $(CFLAGS) $(LIB) collector_pmc.c  $(INCLUDE) -o collector_pmc

repair_pmc: repair_pmc.c common
	$(CC) $(CFLAGS) $(LIB) repair_pmc.c  $(INCLUDE) -o repair_pmc

tests: iotest.c common
	$(CC) $(CFLAGS) $(LIB) iotest.c  $(INCLUDE) -o iotest

common:
	$(CC) $(CFLAGS) -c common.c $(INCLUDE)
	make -C .. # make the core pmc libraries

clean:
	rm -f G.code
	rm -f $(TARGETS) *.o
	rm -rf *.dSYM
