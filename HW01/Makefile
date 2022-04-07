target: Trap.bin SPO.bin MPIOut.bin Simp.bin

Trap.bin: Trap.c
	mpicc -g -Wall -o Trap.bin Trap.c -lm

SPO.bin: SPO.c
	mpicc -g -Wall -o SPO.bin SPO.c

MPIOut.bin: MPIOut.c
	mpicc -g -Wall -o MPIOut.bin MPIOut.c

Simp.bin: Simp.c
	mpicc -g -Wall -o Simp.bin Simp.c -lm

clean:
	rm -rf *.bin *.o

