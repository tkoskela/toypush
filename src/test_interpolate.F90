program test_interpolate

  use params
  use interpolate
  use initmodule
  use particle
  
  implicit none

  integer :: nr,nz,ir,iz,err

  double precision, dimension(4,veclen) :: y
  double precision, dimension(3,veclen) :: efield,bfield
  double precision, dimension(3,3,veclen) :: jacb 

  type(particle_data) :: prt
  
  err = init(prt)
  
  open(unit=11,file='bfield.dat',action='write')
  open(unit=12,file='efield.dat',action='write')
  
  do ir = 1,veclen
     y(1,:) = rmin + dble(ir)/dble(veclen) * (rmax-rmin)
     do iz = 1,veclen
        y(3,iz) = zmin + dble(iz-1)/dble(veclen-1) * (zmax-zmin)
     end do
     err = b_interpol_analytic(y,bfield,jacb)
     err = e_interpol_tri(y,efield)
     
     do iz = 1,veclen
        write(11,*) y(1,iz),y(3,iz),bfield(1,iz),bfield(2,iz),bfield(3,iz)
        write(12,*) y(1,iz),y(3,iz),efield(1,iz),efield(2,iz),efield(3,iz)
     end do
  end do
  
end program test_interpolate
