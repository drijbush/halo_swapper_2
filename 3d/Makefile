.SUFFIXES: 
.SUFFIXES: .f90 .c


PROG =	halo

SRCS =	constants.f90 halo.f90 halo_base_class.f90 halo_parallel_2.f90 \
	swap_module.f90

OBJS =	constants.o halo.o halo_base_class.o halo_parallel_2.o swap_module.o

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

halo.o: swap_module.o halo_parallel_2.o
halo_base_class.o: constants.o
halo_parallel_2.o: constants.o halo_base_class.o swap_module.o
swap_module.o: constants.o 
