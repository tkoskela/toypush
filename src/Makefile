FC = ftn
OBJS = params.o grid.o particle.o particleIO.o equationOfMotion.o interpolate.o rungekutta.o init.o 

FFLAGS = -g -fopenmp

DFLAGS = -DDEBUG

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
