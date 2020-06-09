# Makefile for TSP challenge code

# Do "make depend" to establish the right file dependencies.
# Do "make" or "make all" to compile and link the programs.

# To enable internal timing within the greedy heuristic, add -DTIMING
# as an option on the CFLAGS line.  This is an option because the 
# required system calls are not available on all systems.

OBJS = gmain.o Instance.o Utilities.o
LOBJS = length.o Instance.o Utilities.o

CFILES = gmain.c Instance.c Utilities.c length.c portgen.c portcgen.c \
	portmgen.c

CFLAGS = -O

all:	greedy length portgen portcgen portmgen

greedy:	$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -lm -o greedy

length:	$(LOBJS)
	$(CC) $(CFLAGS) $(LOBJS) -lm -o length

portgen: portgen.o
	$(CC) $(CFLAGS) portgen.o -o portgen

portcgen: portcgen.o
	$(CC) $(CFLAGS) portcgen.o -lm -o portcgen

portmgen: portmgen.o
	$(CC) $(CFLAGS) portmgen.o -o portmgen

depend:
	makedepend $(CFILES)

clean:
	rm -f *.o *.bak greedy length portgen portcgen portmgen *~
