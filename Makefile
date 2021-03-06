PROG =	halo

SRCS =	halo.f90 swap_module.f90

OBJS =	halo.o swap_module.o

LIBS =	

CC = cc
CFLAGS = -O
FC = f77
FFLAGS = -O
F90 = mpif90
F90FLAGS = -O -g -std=f2008 -Wall -Wextra -fcheck=all -finit-real=snan -Wimplicit-interface -Wuse-without-only
LDFLAGS = 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

halo.o: swap_module.o
swap_module.o: 
