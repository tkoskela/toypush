!> This module stores relevant particle data
!<
module particle

  use params
  
  implicit none

  integer, parameter :: ALREADYALLOCATED = params_particleModuleErrorId + 1
  integer, parameter :: NOTALLOCATED     = params_particleModuleErrorId + 2
  
  type particle_data

     double precision, allocatable, dimension(:) :: mass
     double precision, allocatable, dimension(:) :: charge
     double precision, allocatable, dimension(:,:) :: rpz
     double precision, allocatable, dimension(:) :: mu
     double precision, allocatable, dimension(:) :: rho_par

     logical :: isAllocated
     
  end type particle_data

contains

  !> Updates the phase variables rpz and rho_par of the particle
  !<
  function particle_updatePhase(prt,newPhase,veclen,iblock) result(err)

    integer, intent(in) :: veclen !> block size
    integer, intent(in) :: iblock !> block id  
    type(particle_data), intent(inout) :: prt
    double precision, dimension(4,veclen), intent(in) :: newPhase

    integer :: err
    integer :: iv,iglob

    err = 0
    
    do iv = 1,veclen
       iglob = (iblock - 1) * veclen + iv
       prt%rpz(:,iglob)   = newPhase(1:3,iv)
       prt%rho_par(iglob) = newPhase(4,iv)
    end do
    
  end function particle_updatePhase

  !> Returns the phase variables rpz and rho_par of the particle
  !<
  function particle_getPhase(prt,phase,veclen,iblock) result(err)

    integer, intent(in) :: veclen !> block size
    integer, intent(in) :: iblock !> block id  
    type(particle_data), intent(in) :: prt
    double precision, dimension(4,veclen), intent(out) :: phase

    integer :: err
    integer :: iv, iglob

    err = 0
    
    do iv = 1,veclen
       iglob = (iblock - 1) * veclen + iv
       phase(1:3,iv) = prt%rpz(:,iglob)
       phase(4,iv)   = prt%rho_par(iglob)
    end do
    
  end function particle_getPhase
  
  function particle_init(prt,arraydim) result(err)

    integer, intent(in) :: arraydim !> size of arrays
    type(particle_data), intent(inout) :: prt !> particle data struct
    integer :: err

    if(prt%isallocated) then
       err = ALREADYALLOCATED
       return
    end if

    allocate(prt%rpz(3,arraydim))
    allocate(prt%mass(arraydim))
    allocate(prt%charge(arraydim))
    allocate(prt%mu(arraydim))
    allocate(prt%rho_par(arraydim))

    prt%isAllocated = .true.
    err = 0
    
  end function particle_init

  function particle_deallocate(prt) result(err)

    type(particle_data), intent(inout) :: prt
    integer :: err

    if(.not. prt%isallocated) then
       err = NOTALLOCATED
       return
    end if

    deallocate(prt%rpz)
    deallocate(prt%mass)
    deallocate(prt%charge)
    deallocate(prt%mu)
    deallocate(prt%rho_par)

    prt%isallocated = .false.
    err = 0
    
  end function particle_deallocate
  
end module particle
