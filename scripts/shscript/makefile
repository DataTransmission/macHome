# Names and option to fortran compiler
# Compaq compiler and its options
F90 = ifort 
FFLAGS = -g -O5

#NCLIBPATH = -L/nethome/gchen/netcdflib/lib
#NCIC = -I/nethome/gchen/netcdflib/include
#NCLIB     = -lnetcdf
NCLIBPATH = -L/share/apps/netcdf/4.0.1/intel/64/lib
NCIC = -I/share/apps/netcdf/4.0.1/intel/64/include
NCLIB     = -lnetcdf
#FFLAGS = -g2 -O0 -ladebug -free -check bounds -check overflow -fpe

#IBM
#F90 = xlf90 -qfree=f90
#FFLAGS = -O3 -qstrict
#FFLAGS = -g -C -qflttrap=ov:zero:inv:en -qfree=f90

# HP Exemplar compiler and its options
#F90 = f90 +source=free 
#FFLAGS = -g +fp_exception +check=all
#FFLAGS = -O2 
#FFLAGS = -f free -O3,aggress,unroll2,pipeline2 -apad -X 32

# SGI compiler and its options
#F90 = f90 -freeform -mips4 -n32
#FFLAGS = -O3 -mips4
#FFLAGS = -g -C -freeform -64 -mips4 -DEBUG:trap_uninitialized=ON:verbose_runtime=ON:div_check=3

# SUN compiler and its options
#F90 = f90 -free
#FFLAGS = -fast -native -dalign 
#FFLAGS = -g -dalign  -ftrap=%all,no%inexact,no%underflow

# precedence rules
.SUFFIXES: .o .f90
%.o:%.f90
	$(F90)  -c $(FFLAGS) $(NCIC) $< 
#.f90.o:
#	$(F90) -c $(FFLAGS) $<

SOURCES=  myncread.f90 tcdetect.f90
#
# this is the list of object files needed to build the executable
#
tcOBJS= myncread.o tcdetect.o

tcdetect: $(tcOBJS)
	$(F90) $(FFLAGS) $(NCIC) $(SOURCES) -o tcdetect $(NCLIBPATH) $(NCLIB) 
     

depend:
	sfmakedepend $(NCIC) $(SOURCES) 

make clean:
	rm *.o *.mod

# DO NOT DELETE THIS LINE - used by make depend

tcdetect.o: myncread.o

