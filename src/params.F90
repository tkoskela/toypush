!> This module contains global parameters
!<
module params

  implicit none

  ! module ids
  integer, parameter :: params_particleModuleErrorId    = 100
  integer, parameter :: params_eomModuleErrorId         = 200
  integer, parameter :: params_interpolateModuleErrorId = 300

  ! physics parameters
  double precision, parameter :: pi = 3.14159265358979D0
  double precision, parameter :: twopi = 6.28318530717959D0
  double precision, parameter :: speedoflight = 299792458.0D0
  double precision, parameter :: protonmass = 1.6720D-27
  double precision, parameter :: electronmass = 9.1094D-31
  double precision, parameter :: unitcharge = 1.6022D-19
  double precision, parameter :: epsilon0 = 8.8542D-12

  ! run control parameters
  double precision, parameter :: rmin = 1D0
  double precision, parameter :: rmax = 3D0
  double precision, parameter :: zmin = -1D0
  double precision, parameter :: zmax = 1D0 
  integer, parameter :: veclen = 2**6
  integer, parameter :: params_nprtPerRank = 2**15
  integer, parameter :: nt = 1000
#ifdef MULTIPLEELEMENTS
  integer, parameter :: params_nnode = 4
  integer, parameter :: params_ntri = 2
#else
  integer, parameter :: params_nnode = 3
  integer, parameter :: params_ntri = 1
#endif
  double precision, parameter :: dt = 2.0D-8

  integer, save :: params_nprt
  
end module params
