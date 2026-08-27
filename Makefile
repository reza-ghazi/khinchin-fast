SHELL := /bin/bash
CC ?= cc
CFLAGS ?= -O3 -march=native -flto -fopenmp -DNDEBUG -Wall -Wextra -Wpedantic
CPPFLAGS ?=
LDLIBS ?= -fopenmp -lflint -lmpfr -lgmp -lm -lpthread

.PHONY: all clean check bench

all: khinchin-fast

# Verify the C program against the reference digits and both backends,
# then every port with an installed toolchain against the C output.
check: khinchin-fast
	./khinchin-fast 100 /tmp/khinchin-check-100.txt
	@diff <(tr -d '\n' < /tmp/khinchin-check-100.txt) \
	      <(tr -d '\n' < khinchin-100.txt) \
	  && echo "C 100-digit reference: OK"
	./khinchin-fast 1000 /tmp/khinchin-check-scp.txt
	KHINCHIN_BACKEND=arb ./khinchin-fast 1000 /tmp/khinchin-check-arb.txt
	@diff /tmp/khinchin-check-scp.txt /tmp/khinchin-check-arb.txt \
	  && echo "C scp/arb cross-check: OK"
	@rm -f /tmp/khinchin-check-*.txt
	./check-ports.sh

bench: khinchin-fast
	./bench.sh

khinchin-fast: khinchin_fast.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f khinchin-fast
