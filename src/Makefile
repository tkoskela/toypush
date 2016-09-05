FC = gfortran
OBJS = params.o grid.o particle.o particleIO.o equationOfMotion.o interpolate.o rungekutta.o toypush.o

FLAGS = -g -ffpe-trap=zero,invalid,overflow,underflow -fno-backtrace

toypush : $(OBJS)
	$(FC) $(FLAGS) -o $@ $^

%.o : %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o *.mod toypush
