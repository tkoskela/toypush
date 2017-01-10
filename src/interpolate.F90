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

  function e_interpol_tri(y,itri,evec) result(err)

    use grid_module
    
    double precision, intent(in),  dimension(4,veclen) :: y    !> R,phi,z,rho_parallel
    integer, intent(in), dimension(veclen) :: itri
    double precision, intent(out), dimension(3,veclen) :: evec !> Er,Ephi,Ez
    integer :: err

    integer :: iv
    double precision, dimension(3) :: bc_coords
    double precision, dimension(2) :: dx
    integer :: inode, icomp

    err = 0
    evec = 0D0
    do iv = 1,veclen

       dx(1) = y(1,iv) - grid_mapping(1,3,itri)
       dx(2) = y(3,iv) - grid_mapping(2,3,itri)
       bc_coords(1:2) = grid_mapping(1:2,1,itri) * dx(1) + grid_mapping(1:2,2,itri) * dx(2)
       bc_coords(3) = 1.0D0 - bc_coords(1) - bc_coords(2)
            
       do inode = 1,3
          do icomp = 1,3
             evec(icomp,iv) = evec(icomp,iv) + grid_efield(icomp,inode) * bc_coords(inode)
          end do
       end do

       !write(313,*) sngl(y(1,iv)),sngl(y(3,iv)),sngl(evec(:,iv))
       
    end do
    
  end function e_interpol_tri

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
  !> Returns values from an analytic function psi = (1-(2-R)^2)*(1-z^2),
  !> Br = 1/r dpsi/dz, Bz = -1/r dpsi/dR
  !> Input R should be within 1 and 3, input z should be within -1 and 1.
  !<
  function b_interpol_analytic(y,bvec,jacb) result(err)

    double precision, parameter :: br0 = 1D-2
    double precision, parameter :: bz0 = 1D-2
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

    double precision :: over_r, over_r2, r, z
    
    err = 0
    do iv = 1,veclen
       r = y(1,iv) - 2D0
       z = y(3,iv)
       
       over_r = 1.0D0 / y(1,iv)
       over_r2 = over_r / y(1,iv)
       
       bvec(1,iv) = br0 * 2D0 * (r ** 2 - 1D0) * z
       bvec(2,iv) = bp0 * over_r
       bvec(3,iv) = bz0 * 2D0 * (1D0 - z ** 2) * r

       jacb(1,1,iv) = br0 * 4D0 * r * z
       jacb(1,2,iv) = 0D0
       jacb(1,3,iv) = br0 * 2D0 * (r ** 2 - 1D0)

       jacb(2,1,iv) = -1D0 * bp0 * over_r2
       jacb(2,2,iv) = 0D0
       jacb(2,3,iv) = 0D0

       jacb(3,1,iv) = bz0 * 2D0 * (1D0 - z ** 2)
       jacb(3,2,iv) = 0D0
       jacb(3,3,iv) = -bz0 * 4D0 * r * z
    end do

  end function b_interpol_analytic
  
end module interpolate
