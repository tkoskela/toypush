FFLAGS = -g -O3

ifeq ($(NERSC_HOST),edison)
	FC = ftn
else ifeq ($(NERSC_HOST),cori)
	FC = ftn
else ifeq ($(NERSC_HOST),carl)
	FC = mpiifort
	FFLAGS += -xMIC-AVX512
else
	FC = mpifort
endif

OBJS = params.o grid.o particle.o particleIO.o equationOfMotion.o interpolate.o rungekutta.o init.o 
OBJS_AOS = params.o grid.o particle_aos.o particleIO_aos.o equationOfMotion_aos.o interpolate_aos.o rungekutta_aos.o init_aos.o 

DFLAGS += -DMPI
#DFLAGS += -DDEBUG
#DFLAGS += -DOPENMP
#DFLAGS += -DEXPERIMENT

toypush : $(OBJS_AOS) toypush.o
	$(FC) $(FFLAGS) $(DFLAGS) -o $@ $^

test_interpolate : $(OBJS) test_interpolate.o
	$(FC) $(FFLAGS) $(DFLAGS) -o $@ $^

etags :
	etags *.F90

%.o : %.F90
	$(FC) $(FFLAGS) $(DFLAGS) -c $<

clean:
	rm -f *.o *.mod toypush
