!> This program is a simple implementation of a RK4 particle push routine to test theoretical peak performance
!> @author T. Koskela
!> @date Aug 17 2016
!<
program toypush

  use params
  use rk4, only: rk4_push
  use particle, only: particle_data
  use initmodule, only : init

#ifdef OPENMP
  use omp_lib
#endif
#ifdef MPI
  use mpi
#endif
  
  implicit none  
  
  integer :: err
  type(particle_data) :: prt
  integer :: it,iblock
  integer :: nblock
  
  integer :: pid

  integer :: num_procs, my_id, ith

  double precision :: t1,t2

  my_id = 0

#ifdef MPI
  call mpi_init(err)
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, err)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, err)
#endif

#ifdef OPENMP
  !$omp parallel
  !$omp master
  if(my_id .eq. 0) then
     write(*,*) 'number of OpenMP threads = ',omp_get_num_threads()
  end if
  !$omp end master
  !$omp end parallel
#endif

  if(my_id .eq. 0) write(*,*) 'program toypush started'
  if(my_id .eq. 0) write(*,*) 'veclen = ',veclen
  if(my_id .eq. 0) write(*,*)
  
  if(my_id .eq. 0) write(*,*) 'initializing particles with ',params_nprt,'particles.'
  if(my_id .eq. 0) write(*,*) 'initializing grid      with ',params_nnode,'nodes.'
  if(my_id .eq. 0) write(*,*) 'initializing grid      with ',params_ntri,'triangles.'
  err = init(prt)
  if(err .ne. 0) stop
  if(my_id .eq. 0) write(*,*) 'done initialising'
  if(my_id .eq. 0) write(*,*)
  
  nblock = params_nprt / veclen  
  if(my_id .eq. 0) write(*,*) 'pushing particles in ',nblock,' blocks'
  
  !pid = 88
  !open(unit=15,file='orbit.dat',action='write')

#ifdef OPENMP
  t1 = omp_get_wtime()
#else
  call cpu_time(t1)
#endif

  !$omp parallel do private(iblock, it)
  do iblock = 1,nblock

#ifdef MPI
     if (mod(iblock,num_procs) .eq. my_id) then
#endif
     
        do it = 1,nt
           err = rk4_push(prt, iblock)
           
#ifdef VERBOSE
           !$omp critical
           if (mod(it,nt) .eq. 0) then
              ith = omp_get_thread_num()
              if(my_id .eq. 0) write(*,*) 'Thread ',ith,' completed block ',iblock
           end if
           !$omp end critical
#endif
        end do
#ifdef MPI
     end if
#endif
  end do
  !close(15)

#ifdef OPENMP
  t2 = omp_get_wtime()
#else
  call cpu_time(t2)
#endif

  if(my_id .eq. 0) write(*,*) 'done pushing'
  if(my_id .eq. 0) write(*,*) 'spent ',t2-t1,'s'
  if(my_id .eq. 0) write(*,*)

  !call mpi_gatherv
  
  if(my_id .eq. 0)  write(*,*) 'finalizing'
  err = finalize(prt)
#ifdef MPI
  call mpi_finalize(err)
#endif
  if(my_id .eq. 0) write(*,*) 'done finalizing'
  if(my_id .eq. 0) write(*,*)
  
  if(my_id .eq. 0) write(*,*) 'program toypush successfully completed!'
  
contains
  
  function finalize(prt) result(err)

    use rk4, only: rk4_deallocate
    use particle, only : particle_data, particle_deallocate
#ifdef USEIO
    use particleio, only : particleio_write
#endif
    
    implicit none

    type(particle_data) :: prt
    integer :: err

#ifdef USEIO
    err = particleio_write(prt, 'endstate.dat')
#endif
    
    err = particle_deallocate(prt)
    err = rk4_deallocate()
    
  end function finalize
  
end program toypush
