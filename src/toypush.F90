!> This program is a simple implementation of a RK4 particle push routine to test theoretical peak performance
!> @author T. Koskela
!> @date Aug 17 2016
!<
program toypush

  use params
  use rk4, only: rk4_push
  use particle, only: particle_data
  use initmodule, only : init
  
  implicit none  
  
  integer :: err
  type(particle_data) :: prt
  integer :: it

  integer :: pid

  write(*,*) 'program toypush started'
  write(*,*)
  
  write(*,*) 'initializing particles'
  err = init(prt)
  write(*,*) 'done initialising'
  write(*,*)
  
  write(*,*) 'pushing particles'
  
  !pid = 88
  !open(unit=15,file='orbit.dat',action='write')
  
  !$omp parallel do private(it)
  do it = 1,nt
     err = rk4_push(dt, prt)

     !write(15,*) prt%rpz(:,pid)
     if (mod(it,nt/10) .eq. 0) then
        write(*,*) 'completed time step ',it,' out of ',nt
     end if
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
