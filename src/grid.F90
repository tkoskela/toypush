module grid_module

  implicit none

  integer, parameter :: nnode = 3
  integer, parameter :: ntri = 1

  double precision, dimension(3,nnode) :: grid_efield
  double precision, dimension(2,3,nnode) :: grid_mapping
  double precision, private, dimension(2,nnode) :: grid_node
  integer, private, dimension(3,ntri) :: grid_tri

contains
  
  function get_coords(inode,coords) result(err)

    integer, intent(in) :: inode
    double precision, dimension(2), intent(out) :: coords
    integer :: err

    err = 0

    if(inode .gt. nnode) then
       err = 1
       return
    end if
    
    coords = grid_node(:,inode)
    
  end function get_coords

  function get_nodes(itri,nodes) result(err)

    integer, intent(in) :: itri
    integer, dimension(3), intent(out) :: nodes
    integer :: err

    err = 0
    
    if(itri .gt. ntri) then
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
    
    if(itri .gt. ntri) then
       err = 1
       return
    end if

    efield = grid_efield(:,itri)
    
  end function get_efield

  function grid_init_coords(coords) result(err)

    double precision, dimension(2,nnode), intent(in) :: coords
    integer :: err

    double precision :: det
    
    err = 0
    
    grid_tri(:,1) = [1,2,3]
    grid_node = coords

    det = 1D0 / ((coords(2,2)-coords(2,3))*(coords(1,1)-coords(1,3)) &
         + (coords(1,3)-coords(1,2))*(coords(2,1)-coords(2,3)))
    
    grid_mapping(1,1,1) =   det * ( coords(2,2) - coords(2,3) )
    grid_mapping(1,2,1) = - det * ( coords(1,2) - coords(1,3) )
    grid_mapping(2,1,1) = - det * ( coords(2,1) - coords(2,3) )
    grid_mapping(2,2,1) =   det * ( coords(1,1) - coords(1,3) )
    grid_mapping(:,3,1) = coords(:,3)    
    
  end function grid_init_coords

  function grid_init_efield(efield) result(err)

    double precision, dimension(3,nnode), intent(in) :: efield
    integer :: err
    
    err = 0
    
    grid_efield = efield
    
  end function grid_init_efield
  
end module grid_module
