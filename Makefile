# Compiler/linker setup ------------------------------------------------------

PLATFORM = linux
CC       = nvc
CFLAGS   = -acc -fast -gpu=cc120 -Minfo=accel -std=c11
OSLIBS   =
LDFLAGS  = 

# make all: all OBJS are compiled and executable available
# make <individual>
# make clean

OBJS = msa
all: $(OBJS)

msa: msa.c msa_fun.o
	$(CC) $(CFLAGS) msa.c msa_fun.o $(LDFLAGS) -o msa

# library
msa_fun.o: msa_fun.c msa_fun.h
	$(CC) $(CFLAGS) -c msa_fun.c $(LDFLAGS)

# Maintenance and stuff ------------------------------------------------------
clean:
	rm -f $(OBJS) *.o core