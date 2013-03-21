
CC      = gcc
CFLAGS  = -g -Wall -std=c11 -O2 -fopenmp
LDFLAGS = -lrt -lhdf5_hl -lhdf5 -llapacke -lgsl -lgslcblas -lfftw3_omp -lfftw3 -lm

discrete-convolutions.o: discrete-convolutions.c

# test.o: test.c
# 	$(CC) $(CFLAGS) -o test.o -c test.c 

# test: test.o discrete-convolutions.o
# 	$(CC) $(CFLAGS) $(LIBS) -o test test.o discrete-convolutions.o 


test-multivariate-normal: multivariate-normal.o utilities.o
multivariate-normal.o: utilities.o

test-generalised-inverse-gaussian: generalised-inverse-gaussian.o

test-multivariate-generalised-hyperbolic: multivariate-generalised-hyperbolic.o generalised-inverse-gaussian.o multivariate-normal.o utilities.o

test-univariate-generalised-hyperbolic: multivariate-generalised-hyperbolic.o generalised-inverse-gaussian.o multivariate-normal.o utilities.o

test-trawl-process: trawl-process.o multivariate-generalised-hyperbolic.o generalised-inverse-gaussian.o multivariate-normal.o utilities.o

test-argp: test-argp.o

test-utilities: utilities.o

simulate-trawl-process: trawl-process.o multivariate-generalised-hyperbolic.o generalised-inverse-gaussian.o multivariate-normal.o utilities.o

simulate-homogeneous-levy-basis: multivariate-normal.o utilities.o

simulate-vector-field: discrete-convolutions.o utilities.o

README: README.md
	markdown README.md > README.html

install: simulate-trawl-process README
	cp simulate-trawl-process software/
	cp README.md              software/
	cp README.html            software/

clean:
	rm -f *.o *~ 

dist-clean:
	rm -f *.o *~ software/*
