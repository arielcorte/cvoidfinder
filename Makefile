# Makefile for libvoidfinder (pure C, OpenMP optional)
#
#   make            build static library + example
#   make OPENMP=1   build with OpenMP parallelism
#   make example    build only the example program
#   make clean

CC      ?= cc
AR      ?= ar
CSTD    ?= -std=c11
CFLAGS  ?= -O2 -Wall -Wextra $(CSTD) -Iinclude -Isrc
LDFLAGS ?=
LDLIBS  ?= -lm

ifeq ($(OPENMP),1)
  CFLAGS  += -fopenmp
  LDFLAGS += -fopenmp
endif

# getline() needs POSIX on glibc
CFLAGS += -D_POSIX_C_SOURCE=200809L

SRC     := src/rng.c src/catalogue.c src/chainmesh.c src/voidfinder.c
OBJ     := $(SRC:.c=.o)
LIB     := libvoidfinder.a
EXAMPLE := void_finder_example

.PHONY: all example clean

all: $(LIB) example

$(LIB): $(OBJ)
	$(AR) rcs $@ $(OBJ)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

example: $(EXAMPLE)

$(EXAMPLE): examples/void_finder_example.c $(LIB)
	$(CC) $(CFLAGS) $< -o $@ -L. -lvoidfinder $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OBJ) $(LIB) $(EXAMPLE)
