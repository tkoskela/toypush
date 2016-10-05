# toypush
Simple PIC particle push program. Initialises a group of 32 768 particles and pushes them for 1000 steps using a gyrokinetic EOM [Chang et al. PoP 16 056108 2009]. Magnetic field values are evaluated from an analytic function and electric field values are linearly interpolated from three grid points using barycentric coordinates.

# Build Instructions
cd src
make

The Makefile detects the fortran mpi compiler on NERSC platforms Edison, Cori and Carl. On other platforms it tries to use mpifort by default. To use another compiler, edit the Makefile on line 11. Parallelisation by MPI is enabled by default. To use OpenMP, add -DOPENMP to DFLAGS. NOTE: scaling is poor with OpenMP, probably due to false sharing of grid data.

# Run Instructions
Run the executable toypush anywhere. No input files are needed. Run parameters are hard-coded in params.F90.