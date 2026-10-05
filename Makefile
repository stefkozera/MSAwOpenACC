# Compiler/linker setup ------------------------------------------------------

# Mac OS X-specific flags.  Comment these out if using Linux.
#PLATFORM = osx
#CC       = gcc
#CFLAGS   = -fast -Wall
#OSLIBS   = -Wl,-framework -Wl,IOKit
#LDFLAGS  =

# Linux-specific flags.  Comment these out if using Mac OS X.
PLATFORM = linux
CC       = nvc
CFLAGS   = -acc -fast -gpu=cc120 -Minfo=accel -std=c11
OSLIBS   =
LDFLAGS  = 

#-std=c++11: special modifier to allow for lambda expressions in TBB

# make all: all OBJS are compiled and executable available
# make <individual>
# make clean
# '-ltbb' has to be included at the end of the command (can't be in the CFLAGS)
# do not use spaces (tabs are allowed)
# Example programs -----------------------------------------------------------

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