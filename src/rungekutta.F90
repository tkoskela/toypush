module rk4

  use params
  
  implicit none  
  
  double precision, allocatable, dimension(:,:) :: ytmp,dy,dyt,dym,y,y2
  double precision, allocatable, dimension(:,:) :: Efield !> Electric field vector (Er,Ephi,Ez)
  double precision, allocatable, dimension(:,:) :: Bfield !> Magnetic field vector (Br,Bphi,Bz)
  !> Magnetic field jacobian
  !> |  dBRdR,   dBRdphi,   dBRdz  |
  !> | dBphidR, dBphidphi, dBphidz |
  !> |  dBzdR,   dBzdphi,   dBzdz  |
  !<
  double precision, allocatable, dimension(:,:,:) :: jacB
  
contains
  
  function rk4_push(prt, iblock) result(err)

    use eom, only : eom_eval
    use interpolate, only : b_interpol_analytic, e_interpol_tri
    use particle, only : particle_data, particle_getphase, particle_updatephase
    
    implicit none

    type(particle_data), intent(inout) :: prt
    integer, intent(in) :: iblock
    
    double precision :: hdt
    integer :: err

    integer :: iv,iy

    hdt = dt * 0.5D0

    err = particle_getPhase(prt,y,veclen,iblock)

#ifdef DEBUG
    err = check_bounds(y)
    if(err .eq. 1) stop
#endif
    
    ! get derivs with existing E-field
    err = e_interpol_tri(y,efield)
    err = b_interpol_analytic(y,bfield,jacb)
    err = eom_eval(y   ,bfield,jacb,efield,dt,dy ,prt%mu,prt%charge,prt%mass)

    do iy = 1,4
       do iv = 1,veclen
          ytmp(iy,iv) = y(iy,iv) + hdt * dy(iy,iv)
       end do
    end do
#ifdef DEBUG
    err = check_bounds(ytmp)
    if(err .eq. 1) stop
#endif
    err = e_interpol_tri(ytmp,efield)
    err = b_interpol_analytic(ytmp,bfield,jacb)
    err = eom_eval(ytmp,bfield,jacb,efield,dt,dyt,prt%mu,prt%charge,prt%mass)

    do iy = 1,4
       do iv = 1,veclen
          ytmp(iy,iv) = y(iy,iv) + hdt * dyt(iy,iv)
       end do
    end do
#ifdef DEBUG
    err = check_bounds(ytmp)
    if(err .eq. 1) stop
#endif

    err = e_interpol_tri(ytmp,efield)
    err = b_interpol_analytic(ytmp,bfield,jacb)
    err = eom_eval(ytmp,bfield,jacb,efield,dt,dym,prt%mu,prt%charge,prt%mass)
    
    do iy = 1,4
       do iv = 1,veclen
          ytmp(iy,iv) = y(iy,iv) + dt * dym(iy,iv)
       end do
    end do
#ifdef DEBUG
    err = check_bounds(ytmp)
    if(err .eq. 1) stop
#endif
    do iy = 1,4
       do iv = 1,veclen
          dym(iy,iv) = dyt(iy,iv) + dym(iy,iv)
       end do
    end do

    err = e_interpol_tri(ytmp,efield)
    err = b_interpol_analytic(ytmp,bfield,jacb)
    err = eom_eval(ytmp,bfield,jacb,efield,dt,dyt,prt%mu,prt%charge,prt%mass)
    
    ! Obtain new_phase
    do iy = 1,4
       do iv = 1,veclen
          y2(iy,iv) = y(iy,iv) + dt / 6.0D0 * ( dy(iy,iv) + dyt(iy,iv) + 2.0D0*dym(iy,iv) )
       end do
    end do
#ifdef DEBUG
    err = check_bounds(y2)
    if(err .eq. 1) stop
#endif

    err = particle_updatePhase(prt,y2,veclen,iblock)
    
  end function rk4_push

  function rk4_init() result(err)

    integer :: err
    
    allocate(ytmp(4,veclen))
    allocate(y(   4,veclen))
    allocate(y2(  4,veclen))
    allocate(dy(  4,veclen))
    allocate(dyt( 4,veclen))
    allocate(dym( 4,veclen))

    allocate(efield(3,veclen))
    allocate(bfield(3,veclen))
    allocate(jacb(3,3,veclen))

    err = 0
    
  end function rk4_init

  function rk4_deallocate() result(err)

    integer :: err
    
    deallocate(ytmp)
    deallocate(y)
    deallocate(y2)
    deallocate(dy)
    deallocate(dyt)
    deallocate(dym)

    deallocate(efield)
    deallocate(bfield)
    deallocate(jacb)
    
    err = 0

  end function rk4_deallocate

  function check_bounds(y) result(err)

    use params, only : rmin,rmax,zmin,zmax,veclen

    integer :: err
    double precision, intent(in), dimension(4,veclen) :: y
    
    integer :: iv

    do iv = 1,veclen
       if(y(1,iv) .gt. rmax .or. y(1,iv) .lt. rmin &
            .or. y(3,iv) .gt. zmax .or. y(3,iv) .lt. zmin) then
          err = 1
          write(*,*) 'particle ',iv,' is out of bounds!'
          return
       end if       
    end do
    
  end function check_bounds
  
end module rk4
