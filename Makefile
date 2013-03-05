
CC      = gcc
CFLAGS  = -g -Wall -std=c11 -O2 -fopenmp
LDFLAGS = -lrt -llapacke -lgsl -lgslcblas -lfftw3_omp -lfftw3 -lm

# discrete-convolutions.o: discrete-convolutions.c
# 	$(CC) $(CFLAGS) -o discrete-convolutions.o -c discrete-convolutions.c 

# test.o: test.c
# 	$(CC) $(CFLAGS) -o test.o -c test.c 

# test: test.o discrete-convolutions.o
# 	$(CC) $(CFLAGS) $(LIBS) -o test test.o discrete-convolutions.o 

test-generalised-inverse-gaussian: generalised-inverse-gaussian.o

test-multivariate-normal: multivariate-normal.o utilities.o

test-multivariate-generalised-hyperbolic: \
	multivariate-generalised-hyperbolic.o generalised-inverse-gaussian.o \
	multivariate-normal.o utilities.o
