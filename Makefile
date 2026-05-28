CC      = cc
CFLAGS  = -std=c11 -O2 -Wall -Wextra -Wno-unused-parameter
LDFLAGS = -lm

MOL_SRC = main.c chem.c mol.c mutate.c seeds.c
MOL_OBJ = $(MOL_SRC:.c=.o)

ELEM_SRC = disc_elem.c chem.c nuclear.c
ELEM_OBJ = $(ELEM_SRC:.c=.o)

all: exre exre1

exre: $(MOL_OBJ)
	$(CC) $(CFLAGS) -o $@ $(MOL_OBJ) $(LDFLAGS)

exre1: $(ELEM_OBJ)
	$(CC) $(CFLAGS) -o $@ $(ELEM_OBJ) $(LDFLAGS)

%.o: %.c exre.h nuclear.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o exre exre1 exre.log exre_elements.log

.PHONY: all clean
