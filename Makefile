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

msa: msa.o msa_traceback.o msa_grid.o msa_msa.o
	$(CC) $(CFLAGS) msa.o msa_traceback.o msa_grid.o msa_msa.o $(LDFLAGS) -o msa

# library
msa.o: msa.c msa.h
	$(CC) $(CFLAGS) -c msa.c $(LDFLAGS)

msa_traceback.o: msa_traceback.c msa_traceback.h
	$(CC) $(CFLAGS) -c msa_traceback.c $(LDFLAGS)

msa_grid.o: msa_grid.c msa_grid.h
	$(CC) $(CFLAGS) -c msa_grid.c $(LDFLAGS)

msa_msa.o: msa_msa.c msa_msa.h
	$(CC) $(CFLAGS) -c msa_msa.c $(LDFLAGS)

# Maintenance and stuff ------------------------------------------------------
clean:
	rm -f $(OBJS) *.o core