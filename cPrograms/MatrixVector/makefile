EXECS=MatrixVector
MPICC?=mpicc

all: ${EXECS}

${EXECS}: ${EXECS}.c
	${MPICC} -o ${EXECS} ${EXECS}.c

clean:
	rm -f ${EXECS}
