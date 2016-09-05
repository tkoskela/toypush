!> This program is a simple implementation of a RK4 particle push routine to test theoretical peak performance
!> @author T. Koskela
!> @date Aug 17 2016
!<
program toypush

  use params
  use rk4, only: rk4_push
  use particle, only: particle_data
  use particleIO, only: particleio_write
  use grid_module, only : get_coords, get_nodes, grid_init_coords, grid_init_efield
  
  implicit none  
  
  integer :: err
  type(particle_data) :: prt
  integer :: it

  write(*,*) 'program toypush started'
  write(*,*)
  
  write(*,*) 'initializing particles'
  err = init(prt)
  write(*,*) 'done initialising'
  write(*,*)
  
  write(*,*) 'pushing particles'
  do it = 1,nt
     err = rk4_push(dt, prt)
     write(*,*) 'completed time step ',it,' out of ',nt
  end do
  write(*,*) 'done pushing'
  write(*,*)
  
  write(*,*) 'finalizing'
  err = finalize(prt)
  write(*,*) 'done finalizing'
  write(*,*)
  
  write(*,*) 'program toypush successfully completed!'
  
contains

  function init(prt) result(err)

    use rk4, only: rk4_init
    use particle, only : particle_data, particle_init, particle_updatephase
    
    implicit none

    type(particle_data) :: prt
    integer :: i,ir,ith,nr,nth
    integer :: err

    double precision, dimension(4,nprt) :: y

    double precision, dimension(2,3) :: coords
    double precision, dimension(3,3) :: efield

    double precision :: v
    
    nth = 32
    nr = nprt / nth
    
    err = particle_init(prt,nprt)

    v = 1D4
    
    do ir = 1,nr
       do ith = 1,nth
          i = (ith - 1) * nr + ir
          y(1,i) = 1D0 + cos(dble(ith-1)/dble(nth-1) * twopi) * dble(ir)/dble(nr) * 0.8D0
          y(2,i) = 0.0D0
          y(3,i) = sin(dble(ith-1)/dble(nth-1) * twopi) * dble(ir)/dble(nr) * 0.8D0
          y(4,i) = protonmass / unitcharge * v
       end do
    end do

    err = particle_updatePhase(prt,y,nprt)

    prt%mass = 2.0D0 * protonmass
    prt%charge = unitcharge
    prt%mu = 0.5D0 * protonmass * v ** 2

    err = rk4_init()

    coords(:,1) = [rmin,       zmin      ]
    coords(:,2) = [2D0 * rmax, zmin      ]
    coords(:,3) = [rmin,       2D0 * zmax]
    
    err = grid_init_coords(coords)

    !              ER   Ep   Ez
    efield(:,1) = [1D0, 0D0, 0D0] ! node 1
    efield(:,2) = [2D0, 0D0, 0D0] ! node 2
    efield(:,3) = [0D0, 0D0, 0D0] ! node 3

    efield = efield * 1D3
    
    err = grid_init_efield(efield)
    
    err = particleio_write(prt,'inistate.dat')
    
  end function init
  
  function finalize(prt) result(err)

    use rk4, only: rk4_deallocate
    use particle, only : particle_data, particle_deallocate
    
    implicit none

    type(particle_data) :: prt
    integer :: err

    err = particleio_write(prt, 'endstate.dat')
    
    err = particle_deallocate(prt)
    err = rk4_deallocate()
    
  end function finalize
  
end program toypush
