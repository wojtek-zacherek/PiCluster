SRCS = $(wildcard *.c)

PROGS = $(patsubst %.c,%,$(SRCS))

all: $(PROGS)

%: %.c
	mpicc $(CFLAGS)  -o $@.bin $< -lm

clean:
	rm -rf *.bin
