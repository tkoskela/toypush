FC = gfortran
OBJS = params.o particle.o equationOfMotion.o interpolate.o rungekutta.o toypush.o

toypush : $(OBJS)
	$(FC) -o $@ $^

%.o : %.F90
	$(FC) -c $<

clean:
	rm -f *.o *.mod toypush
