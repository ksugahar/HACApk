SYSTEM = INTEL

#intel
ifeq ($(SYSTEM),INTEL)
##OPTFLAGS = -O3 -traceback -ip -heap-arrays -qopenmp
OPTFLAGS = -qopenmp -O3 -lmpi -heap-arrays
CC=mpiicc
F90=mpiifort
CCFLAGS = $(OPTFLAGS)
#F90FLAGS = $(OPTFLAGS) -fpp -assume nounderscore -names uppercase
F90FLAGS = $(OPTFLAGS) -fpp
##F90FLAGS = $(OPTFLAGS) -fpp -check all
#F90FLAGS = -fpe0 -traceback -g -CB -assume nounderscore -names lowercase -fpp -check all
#LDFLAGS = -mkl -trace
LDFLAGS = -mkl
endif
INCLUDE += /home/y42000/anaconda3/pkgs/mpich-3.3.2-hc856adb_0/include/mpif.h

LINK=$(F90)

OBJS= HACApK_lib.o m_HACApK_calc_entry_ij.o \
	 m_HACApK_base.o main.o \


TARGET=Biot_Savart.out

.SUFFIXES:
.SUFFIXES: .o .c .f90

all: $(TARGET)

$(TARGET): $(OBJS)
			$(LINK) -o $@ $(OBJS) $(LDFLAGS)

.c.o: *.c
			$(CC) -c $(CCFLAGS) $<
.f90.o: *.f90
#			echo 'f90 complile'
			$(F90) -c $< $(F90FLAGS)
clean:
	rm -f *.o *.mod $(TARGET)

rmod:
	rm -f m_*.o *.mod

