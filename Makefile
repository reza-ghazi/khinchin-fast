CC ?= cc
CFLAGS ?= -O3 -march=native -flto -fopenmp -DNDEBUG -Wall -Wextra -Wpedantic
CPPFLAGS ?=
LDLIBS ?= -fopenmp -lflint -lmpfr -lgmp -lm -lpthread

.PHONY: all clean

all: khinchin-fast

khinchin-fast: khinchin_fast.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f khinchin-fast
