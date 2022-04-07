target: TrapTimed.bin TrapTimedReduced.bin VectorMultScalarV1.bin VectorMultScalarV2.bin VectorMultScalarV3.bin TTRT.bin

TTRT.bin: TTRT.c
	mpicc -g -Wall -o TTRT.bin TTRT.c -lm

TrapTimed.bin: TrapTimed.c
	mpicc -g -Wall -o TrapTimed.bin TrapTimed.c -lm

TrapTimedReduced.bin: TrapTimedReduced.c
	mpicc -g -Wall -o TrapTimedReduced.bin TrapTimedReduced.c -lm

VectorMultScalarV1.bin: VectorMultScalarV1.c
	mpicc -g -Wall -o VectorMultScalarV1.bin VectorMultScalarV1.c -lm

VectorMultScalarV2.bin: VectorMultScalarV2.c
	mpicc -g -Wall -o VectorMultScalarV2.bin VectorMultScalarV2.c -lm

VectorMultScalarV3.bin: VectorMultScalarV3.c
	mpicc -g -Wall -o VectorMultScalarV3.bin VectorMultScalarV3.c -lm

MPIOut.bin: MPIOut.c
	mpicc -g -Wall -o MPIOut.bin MPIOut.c

Simp.bin: Simp.c
	mpicc -g -Wall -o Simp.bin Simp.c -lm

clean:
	rm -rf *.bin *.o

