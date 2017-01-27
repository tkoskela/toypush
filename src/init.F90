module initmodule

  use params
  
  implicit none

contains
  
  !> Initialize the particles, grid, E-field, etc. This stuff could/should be read from an input
  !> file in a future release.
  !<
  function init(prt) result(err)

    use rk4, only: rk4_init
    use particle, only : particle_data, particle_init, particle_updatephase
#ifdef USEIO
    use particleIO, only: particleio_write
#endif
    use grid_module, only : grid_init,get_coords, get_nodes, grid_init_coords, grid_init_efield

    implicit none

    type(particle_data) :: prt
    integer :: i,ir,ith,nr,nth
    integer :: err

    double precision, dimension(4,params_nprt) :: y

    double precision, dimension(2,3) :: coords
    double precision, dimension(3,3) :: efield
    integer, dimension(3,params_ntri) :: tri

    double precision :: v
    
    nth = 32
    nr = params_nprt / nth
    
    err = particle_init(prt,params_nprt)

    v = 1D4
    
    do ir = 1,nr
       do ith = 1,nth
          i = (ith - 1) * nr + ir
          y(1,i) = rmin + (rmax-rmin) / 2D0 + &
               cos(dble(ith-1)/dble(nth-1) * twopi) * dble(ir)/dble(nr) * 0.5D0
          y(2,i) = 0.0D0
          y(3,i) = sin(dble(ith-1)/dble(nth-1) * twopi) * dble(ir)/dble(nr) * 0.5D0
          y(4,i) = protonmass / unitcharge * v
       end do
    end do

    err = particle_updatePhase(prt,y,params_nprt,1)

    prt%mass = 2.0D0 * protonmass
    prt%charge = unitcharge
    prt%mu = 0.5D0 * protonmass * v ** 2
    
    err = grid_init(params_nnode,params_ntri)

    err = rk4_init()

    coords(:,1) = [rmin,       zmin      ]
    coords(:,2) = [2D0 * rmax, zmin      ]
    coords(:,3) = [rmin,       2D0 * zmax]
    
    tri(:,1) = [1,2,3]

    err = grid_init_coords(coords,tri)

    !              ER   Ep   Ez
    efield(:,1) = [ 0D0, 0D0, 0D0] ! node 1
    efield(:,2) = [ 1D0, 0D0, 0D0] ! node 2
    efield(:,3) = [-1D0, 0D0, 0D0] ! node 3

    efield = efield * 1D0
    
    err = grid_init_efield(efield)
    
#ifdef USEIO
    err = particleio_write(prt,'inistate.dat')
#endif
    
  end function init
end module initmodule
