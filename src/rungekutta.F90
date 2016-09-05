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
  
  function rk4_push(dt,prt) result(err)

    use eom, only : eom_eval
    use interpolate, only : b_interpol_analytic, e_interpol_tri
    use particle, only : particle_data, particle_getphase, particle_updatephase
    
    implicit none

    type(particle_data), intent(inout) :: prt
    double precision, intent(in) :: dt
    
    double precision :: hdt
    integer :: err

    hdt = dt * 0.5D0

    err = particle_getPhase(prt,y,veclen)

#ifdef DEBUG
    err = check_bounds(y)
    if(err .eq. 1) stop
#endif
    
    ! get derivs with existing E-field
    err = e_interpol_tri(y,efield)
    err = b_interpol_analytic(y,bfield,jacb)
    err = eom_eval(y   ,bfield,jacb,efield,dt,dy ,prt%mu,prt%charge,prt%mass)

    ytmp = y + hdt * dy
#ifdef DEBUG
    err = check_bounds(ytmp)
    if(err .eq. 1) stop
#endif
    err = e_interpol_tri(ytmp,efield)
    err = b_interpol_analytic(ytmp,bfield,jacb)
    err = eom_eval(ytmp,bfield,jacb,efield,dt,dyt,prt%mu,prt%charge,prt%mass)

    ytmp = y + hdt * dyt
#ifdef DEBUG
    err = check_bounds(ytmp)
    if(err .eq. 1) stop
#endif

    err = e_interpol_tri(ytmp,efield)
    err = b_interpol_analytic(ytmp,bfield,jacb)
    err = eom_eval(ytmp,bfield,jacb,efield,dt,dym,prt%mu,prt%charge,prt%mass)
    
    ytmp = y + dt * dym
#ifdef DEBUG
    err = check_bounds(ytmp)
    if(err .eq. 1) stop
#endif
    dym = dyt + dym

    err = e_interpol_tri(ytmp,efield)
    err = b_interpol_analytic(ytmp,bfield,jacb)
    err = eom_eval(ytmp,bfield,jacb,efield,dt,dyt,prt%mu,prt%charge,prt%mass)
    
    ! Obtain new_phase
    y2 = y + dt / 6.0D0 * ( dy + dyt + 2.0D0*dym )
#ifdef DEBUG
    err = check_bounds(y2)
    if(err .eq. 1) stop
#endif

    err = particle_updatePhase(prt,y2,veclen)
    
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
