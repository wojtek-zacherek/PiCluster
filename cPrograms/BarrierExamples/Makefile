SRCS = $(wildcard *.c)

PROGS = $(patsubst %.c,%,$(SRCS))

all: $(PROGS)

%: %.c
	mpicc $(CFLAGS)  -o $@.bin $<

clean:
	rm -rf *.bin
