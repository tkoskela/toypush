module particleIO

  use particle
  
  implicit none

  !character(len=12), parameter :: inistatefilename = 'inistate.dat'
  !character(len=12), parameter :: endstatefilename = 'endstate.dat'

  integer, parameter :: chn = 313
  
contains

  !> write out important fields of the particle struct into file
  function particleio_write(prt,fn,ibeg,iend) result(err)

    character(len=*) :: fn
    type(particle_data) :: prt
    integer, intent(in), optional :: ibeg
    integer, intent(in), optional :: iend

    integer :: err,i,i1,i2

    err = 0
    
    if(present(ibeg) .and. present(iend)) then
       i1 = ibeg
       i2 = iend
    else
       i1 = 1
       i2 = size(prt%mass)
    end if
    
    open(unit=chn,file=fn,action='write', iostat=err)
    
    do i = i1,i2
       write(chn,'(5E16.8)') prt % rpz(:,i), prt % mu(i), prt % rho_par(i)
    end do

    close(chn)

  end function particleio_write
  
end module particleIO
