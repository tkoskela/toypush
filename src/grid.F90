module grid_module

  implicit none

  integer :: grid_nnode
  integer :: grid_ntri

  double precision, allocatable, dimension(:,:) :: grid_efield
  double precision, allocatable, dimension(:,:,:) :: grid_mapping
  double precision, allocatable, dimension(:,:) :: grid_node
  integer, allocatable, dimension(:,:) :: grid_tri

contains

  !> Sets grid_nnode and grid_ntri, allocates arrays to the correct size.
  !> Does NOT set values of coordinates or connectivities.
  !<
  function grid_init(nnode,ntri) result(err)

    integer, intent(in) :: nnode,ntri
    integer :: err

    err = 0
    
    if(nnode .ge. 3) then
       grid_nnode = nnode
    else
       err = 1
       return
    end if

    if(ntri .le. nnode .and. ntri .ge. 1) then
       grid_ntri = ntri
    else
       err = 2
       return
    end if

    if(.not. allocated(grid_efield)) then
       allocate(grid_efield(3,nnode))
    else
       err = 3
       return
    end if

    if(.not. allocated(grid_mapping)) then
       allocate(grid_mapping(2,3,nnode))
    else
       err = 4
       return
    end if

    if(.not. allocated(grid_node)) then
       allocate(grid_node(2,nnode))
    else
       err = 5
       return
    end if

    if(.not. allocated(grid_tri)) then
       allocate(grid_tri(3,ntri))
    else
       err = 6
       return
    end if
    
  end function grid_init 
  
  function get_coords(inode,coords) result(err)

    integer, intent(in) :: inode
    double precision, dimension(2), intent(out) :: coords
    integer :: err

    err = 0

    if(inode .gt. grid_nnode) then
       err = 1
       return
    end if
    
    coords = grid_node(:,inode)
    
  end function get_coords

  function get_bc_coords(xy,itri,bc_coords) result(err)

    integer, intent(in) :: itri
    double precision, dimension(2), intent(in) :: xy
    double precision, dimension(3), intent(out) :: bc_coords

    integer :: err
    double precision, dimension(2) :: dx
    
    err = 0
    
    dx(1) = xy(1) - grid_mapping(1,3,itri)
    dx(2) = xy(2) - grid_mapping(2,3,itri)
    bc_coords(1:2) = grid_mapping(1:2,1,itri) * dx(1) + grid_mapping(1:2,2,itri) * dx(2)
    bc_coords(3) = 1.0D0 - bc_coords(1) - bc_coords(2)
    
  end function get_bc_coords

  function get_nodes(itri,nodes) result(err)

    integer, intent(in) :: itri
    integer, dimension(3), intent(out) :: nodes
    integer :: err

    err = 0
    
    if(itri .gt. grid_ntri) then
       err = 1
       return
    end if

    nodes = grid_tri(:,itri)
    
  end function get_nodes

  function get_efield(itri,efield) result(err)

    integer, intent(in) :: itri
    integer, dimension(3), intent(out) :: efield
    integer :: err

    err = 0
    
    if(itri .gt. grid_ntri) then
       err = 1
       return
    end if

    efield = grid_efield(:,itri)
    
  end function get_efield

  function grid_init_coords(coords,tri) result(err)

    double precision, dimension(2,grid_nnode), intent(in) :: coords
    integer, dimension(3,grid_ntri), intent(in) :: tri
    integer :: err

    integer :: inode1,inode2,inode3, itri    
    double precision :: det
    
    err = 0
    
    grid_tri = tri
    grid_node = coords

    do itri = 1,grid_ntri

       inode1 = tri(1,itri)
       inode2 = tri(2,itri)
       inode3 = tri(3,itri)
       
       det = 1D0 / ((coords(2,inode2)-coords(2,inode3))*(coords(1,inode1)-coords(1,inode3)) &
            + (coords(1,inode3)-coords(1,inode2))*(coords(2,inode1)-coords(2,inode3)))
       
       grid_mapping(1,1,itri) =   det * ( coords(2,inode2) - coords(2,inode3) )
       grid_mapping(1,2,itri) = - det * ( coords(1,inode2) - coords(1,inode3) )
       grid_mapping(2,1,itri) = - det * ( coords(2,inode1) - coords(2,inode3) )
       grid_mapping(2,2,itri) =   det * ( coords(1,inode1) - coords(1,inode3) )
       grid_mapping(:,3,itri) = coords(:,inode3)
    end do
    
  end function grid_init_coords

  function grid_init_efield(efield) result(err)

    double precision, dimension(3,grid_nnode), intent(in) :: efield
    integer :: err
    
    err = 0
    
    grid_efield = efield
    
  end function grid_init_efield
  
end module grid_module
