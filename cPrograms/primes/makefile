EXECS=primes
MPICC?=mpicc

all: ${EXECS}

primes: primes.c
	${MPICC} -o primes primes.c

clean:
	rm -f ${EXECS}
