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

  !> Interpolation function for 3-component vector (electric) field from unstructured mesh.
  !<
  function e_interpol_tri(y,itri,evec) result(err)

    use grid_module

    double precision, intent(in),  dimension(veclen,4) :: y    !> R,phi,z,rho_parallel
    integer, intent(in), dimension(veclen) :: itri             !> Grid element index
    double precision, intent(out), dimension(veclen,3) :: evec !> Er,Ephi,Ez          
    integer :: err

    integer :: iv !> Index for vector loop
    double precision, dimension(3) :: bc_coords !> Weight factor for each node
    double precision, dimension(2) :: dx 
    integer, dimension(3) :: nodes !> Node indices for element itri
    integer :: inode, icomp

    integer :: itri_scalar

    err = 0
    evec = 0D0
    if(all(itri .eq. itri(1))) then
       itri_scalar = itri(1)
       !dir$ simd
       !!dir$ nounroll
       !dir$ vector aligned
       do iv = 1,veclen
          
          dx(1) = y(iv,1) - grid_mapping(1,3,itri_scalar)
          dx(2) = y(iv,3) - grid_mapping(2,3,itri_scalar)
          bc_coords(1:2) = grid_mapping(1:2,1,itri_scalar) * dx(1) + grid_mapping(1:2,2,itri_scalar) * dx(2)
          bc_coords(3) = 1.0D0 - bc_coords(1) - bc_coords(2)
          
          do inode = 1,3
             do icomp = 1,3
                evec(iv,icomp) = evec(iv,icomp) + grid_efield(icomp,grid_tri(inode,itri_scalar)) * bc_coords(inode)
             end do
          end do
          
          !write(313,*) sngl(y(1,iv)),sngl(y(3,iv)),sngl(evec(:,iv))
          
       end do
    else
       do iv = 1,veclen
          
          dx(1) = y(iv,1) - grid_mapping(1,3,itri(iv))
          dx(2) = y(iv,3) - grid_mapping(2,3,itri(iv))
          bc_coords(1:2) = grid_mapping(1:2,1,itri(iv)) * dx(1) + grid_mapping(1:2,2,itri(iv)) * dx(2)
          bc_coords(3) = 1.0D0 - bc_coords(1) - bc_coords(2)
          
          nodes = grid_tri(:,itri(iv))
          
          do inode = 1,3
             do icomp = 1,3
                evec(iv,icomp) = evec(iv,icomp) + grid_efield(icomp,nodes(inode)) * bc_coords(inode)
             end do
          end do
          
          !write(313,*) sngl(y(1,iv)),sngl(y(3,iv)),sngl(evec(:,iv))
          
       end do
    end if
    
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
       bvec(iv,:) = 0.0D0
       jacb(iv,:,:) = 0.0D0
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
    
    double precision, intent(in),  dimension(veclen,4) :: y    !> R,phi,z,rho_parallel
    double precision, intent(out), dimension(veclen,3) :: bvec !> Br,Bphi,Bz
    !> Magnetic field jacobian
    !> |  dBRdR,   dBRdphi,   dBRdz  |
    !> | dBphidR, dBphidphi, dBphidz |
    !> |  dBzdR,   dBzdphi,   dBzdz  |
    !<
    double precision, intent(out), dimension(veclen,3,3) :: jacb !> Br,Bphi,Bz          
    integer :: err
    integer :: iv

    double precision :: over_r, over_r2, r, z
    
    err = 0
    do iv = 1,veclen
       r = y(iv,1) - 2D0
       z = y(iv,3)
       
       over_r = 1.0D0 / y(iv,1)
       over_r2 = over_r / y(iv,1)
       
       bvec(iv,1) = br0 * 2D0 * (r ** 2 - 1D0) * z
       bvec(iv,2) = bp0 * over_r
       bvec(iv,3) = bz0 * 2D0 * (1D0 - z ** 2) * r

       jacb(iv,1,1) = br0 * 4D0 * r * z
       jacb(iv,1,2) = 0D0
       jacb(iv,1,3) = br0 * 2D0 * (r ** 2 - 1D0)

       jacb(iv,2,1) = -1D0 * bp0 * over_r2
       jacb(iv,2,2) = 0D0
       jacb(iv,2,3) = 0D0

       jacb(iv,3,1) = bz0 * 2D0 * (1D0 - z ** 2)
       jacb(iv,3,2) = 0D0
       jacb(iv,3,3) = -bz0 * 4D0 * r * z
    end do

  end function b_interpol_analytic
  
end module interpolate
