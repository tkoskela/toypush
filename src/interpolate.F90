!> Module that contains interpolation functions for electric and magnetic fields
!<
module interpolate

  use params
  
  implicit none

contains

  function init_interpol() result(err)

    integer :: err

    err = 0

  end function init_interpol
  
  !> Dummy function that interpolates the electric field vector at a given set of (r,phi,z) coordinates.
  !> Always returns 0
  !<
  function e_interpol_0(y,evec) result(err)

    double precision, intent(in),  dimension(4,veclen) :: y    !> R,phi,z,rho_parallel
    double precision, intent(out), dimension(3,veclen) :: evec !> Er,Ephi,Ez          
    integer :: err

    integer :: iv
    
    err = 0
    do iv = 1,veclen
       evec(:,iv) = 0.0D0
    end do
    
  end function e_interpol_0

  !> Dummy function that interpolates the magnetic field vector at a given set of (r,phi,z) coordinates
  !> Always returns 0
  !<
  function b_interpol_0(y,bvec,jacb) result(err)

    double precision, intent(in),  dimension(4,veclen) :: y    !> R,phi,z,rho_parallel
    double precision, intent(out), dimension(3,veclen) :: bvec !> Br,Bphi,Bz
    !> Magnetic field jacobian
    !> |  dBRdR,   dBRdphi,   dBRdz  |
    !> | dBphidR, dBphidphi, dBphidz |
    !> |  dBzdR,   dBzdphi,   dBzdz  |
    !<
    double precision, intent(out), dimension(3,3,veclen) :: jacb !> Br,Bphi,Bz          
    integer :: err
    integer :: iv

    err = 0
    do iv = 1,veclen
       bvec(:,iv) = 0.0D0
       jacb(:,:,iv) = 0.0D0
    end do
    
  end function b_interpol_0

  !> Function that interpolates the magnetic field vector at a given set of (r,phi,z) coordinates
  !> Returns values from an analytic function psi = (2-R)^2*(1-z)^2, Br = 1/r dpsi/dz, Bz = -1/r dpsi/dR
  !> Input R should be within 0 and 2, input z should be within -1 and 1.
  !<
  function b_interpol_analytic(y,bvec,jacb) result(err)

    double precision, parameter :: br0 = 1D-1
    double precision, parameter :: bz0 = 1D-1
    double precision, parameter :: bp0 = 1D0
    
    double precision, intent(in),  dimension(4,veclen) :: y    !> R,phi,z,rho_parallel
    double precision, intent(out), dimension(3,veclen) :: bvec !> Br,Bphi,Bz
    !> Magnetic field jacobian
    !> |  dBRdR,   dBRdphi,   dBRdz  |
    !> | dBphidR, dBphidphi, dBphidz |
    !> |  dBzdR,   dBzdphi,   dBzdz  |
    !<
    double precision, intent(out), dimension(3,3,veclen) :: jacb !> Br,Bphi,Bz          
    integer :: err
    integer :: iv

    double precision :: over_r, over_r2
    
    err = 0
    do iv = 1,veclen
       over_r = 1.0D0 / y(1,iv)
       over_r2 = over_r / y(1,iv)
       
       bvec(1,iv) = br0 * over_r * 2.0D0 * y(3,iv) * (1.0D0 - y(3,iv)) * (1 - y(1,iv)) ** 2
       bvec(2,iv) = bp0 * over_r
       bvec(3,iv) = bz0 * 2.0D0 * (y(1,iv) - 2.0D0) * (1 - y(3,iv)) ** 2 ! r and 1/r cancel out

       jacb(1,1,iv) = br0 * over_r2 * (y(1,iv) - 2D0) * (y(1,iv) + 2D0) * (y(3,iv) - 1D0) * y(3,iv)
       jacb(1,2,iv) = 0D0
       jacb(1,3,iv) = br0 * over_r  * (2D0 - y(1,iv)**2) * (2D0 * y(3,iv) - 1)

       jacb(2,1,iv) = -1D0 * bp0 * over_r2
       jacb(2,2,iv) = 0D0
       jacb(2,3,iv) = 0D0

       jacb(3,1,iv) = bz0 * (1D0 - y(3,iv)) ** 2
       jacb(3,2,iv) = 0D0
       jacb(3,3,iv) = bz0 * 2D0 * (y(1,iv) - 2D0) * (y(3,iv) - 1D0)
    end do

  end function b_interpol_analytic
  
end module interpolate
