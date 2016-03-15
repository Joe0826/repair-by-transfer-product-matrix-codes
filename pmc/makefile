CC = gcc
CFLAGS = -std=c99 -w -O3 #-g #-DDEBUG #to print debugging info.

LIB = -lgf_complete -lJerasure

all: matrix.c dag_helper.c product_matrix.c MBR_MSR.c rbt.c
	$(CC) $(CFLAGS) -c matrix.c dag_helper.c product_matrix.c MBR_MSR.c rbt.c -I../jerasure-kmg/include

clean:
	rm -f *.o
