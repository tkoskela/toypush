FFLAGS = -g -O3 -traceback -qopenmp

ifeq ($(NERSC_HOST),edison)
	FC = ftn
else ifeq ($(NERSC_HOST),cori)
	FC = ftn
	FFLAGS += -xMIC-AVX512 -qopt-report=5 -qopt-report-phase=vec
else ifeq ($(NERSC_HOST),carl)
	FC = mpiifort
	FFLAGS += -xMIC-AVX512 -qopt-report=5 -qopt-report-phase=vec
else
	FC = mpifort
endif

OBJS = params.o grid.o search.o particle.o particleIO.o equationOfMotion.o interpolate.o rungekutta.o init.o 

#DFLAGS += -DVERBOSE
DFLAGS += -DMPI
DFLAGS += -DMULTIPLEELEMENTS
#DFLAGS += -DDEBUG
#DFLAGS += -DOPENMP
#DFLAGS += -DNODIV
#DFLAGS += -DNOSQRT
#DFLAGS += -DUSEIO

toypush : $(OBJS) toypush.o
	$(FC) $(FFLAGS) $(DFLAGS) -o $@ $^

test_interpolate : $(OBJS) test_interpolate.o
	$(FC) $(FFLAGS) $(DFLAGS) -o $@ $^

etags :
	etags *.F90

%.o : %.F90
	$(FC) $(FFLAGS) $(DFLAGS) -c $<

clean:
	rm -f *.o *.mod toypush
