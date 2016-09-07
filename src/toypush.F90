!> This program is a simple implementation of a RK4 particle push routine to test theoretical peak performance
!> @author T. Koskela
!> @date Aug 17 2016
!<
program toypush

  use params
  use rk4, only: rk4_push
  use particle, only: particle_data
  use initmodule, only : init
  use omp_lib
  
  implicit none  
  
  integer :: err
  type(particle_data) :: prt
  integer :: it,iblock
  integer :: nblock
  
  integer :: pid

  write(*,*) 'program toypush started'
  write(*,*) 'veclen = ',veclen
  write(*,*)

  !$omp parallel
  !$omp master
  write(*,*) 'number of OpenMP threads = ',omp_get_num_threads()
  !$omp end master
  !$omp end parallel
  
  write(*,*) 'initializing',nprt,'particles'  
  err = init(prt)
  write(*,*) 'done initialising'
  write(*,*)
  
  nblock = nprt / veclen  
  write(*,*) 'pushing particles in ',nblock,' blocks'
  
  !pid = 88
  !open(unit=15,file='orbit.dat',action='write')

  !$omp parallel do private(iblock, it)
  do iblock = 1,nblock
     
     do it = 1,nt
        err = rk4_push(prt, iblock)
        
        !write(15,*) prt%rpz(:,pid)
        !$omp critical
        if (mod(it,nt/10) .eq. 0) then
           write(*,*) 'completed time step ',it,' out of ',nt,' in block ',iblock
        end if
        !$omp end critical
     end do
  end do
  !close(15)
  write(*,*) 'done pushing'
  write(*,*)
  
  write(*,*) 'finalizing'
  err = finalize(prt)
  write(*,*) 'done finalizing'
  write(*,*)
  
  write(*,*) 'program toypush successfully completed!'
  
contains
  
  function finalize(prt) result(err)

    use rk4, only: rk4_deallocate
    use particle, only : particle_data, particle_deallocate
    use particleio, only : particleio_write
    
    implicit none

    type(particle_data) :: prt
    integer :: err

    err = particleio_write(prt, 'endstate.dat')
    
    err = particle_deallocate(prt)
    err = rk4_deallocate()
    
  end function finalize
  
end program toypush
