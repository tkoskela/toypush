module search_module

  use params
  
  implicit none

  double precision, parameter :: eps = 1.0D-8

contains

  function search_tr(xy,id) result(err)

    use grid_module, only : grid_mapping, grid_ntri

    double precision, dimension(2) :: xy
    integer :: id

    integer :: itri
    double precision, dimension(3) :: bc_coords !> Weight factor for each node
    double precision, dimension(2) :: dx

    integer :: err

    err = 0
    
    do itri = 1,grid_ntri

       dx(1) = xy(1) - grid_mapping(1,3,itri)
       dx(2) = xy(2) - grid_mapping(2,3,itri)
       bc_coords(1:2) = grid_mapping(1:2,1,itri) * dx(1) + grid_mapping(1:2,2,itri) * dx(2)
       bc_coords(3) = 1.0D0 - bc_coords(1) - bc_coords(2)

       if(minval(bc_coords) .ge. -eps) then
          id = itri
          return
       end if
       
    end do
    
    return
    
  end function search_tr

  function search_tr_vec(y,id) result(err)

    use grid_module, only : grid_mapping, grid_ntri

    double precision, dimension(veclen,4) :: y
    integer, dimension(veclen) :: id

    integer :: itri, iv
    double precision, dimension(3) :: bc_coords !> Weight factor for each node
    double precision :: dx1,dx2 
    double precision, dimension(2) :: dx

    integer :: err
    
    err = 0
    id = -999

    !!$omp simd private(dx,bc_coords)
    !$omp simd 
    !dir$ vector aligned
    vecloop: do iv = 1,veclen
       triangleloop: do itri = 1,grid_ntri
          
          dx(1) = y(iv,1) - grid_mapping(1,3,itri)
          dx(2) = y(iv,3) - grid_mapping(2,3,itri)
          bc_coords(1:2) = grid_mapping(1:2,1,itri) * dx(1) + grid_mapping(1:2,2,itri) * dx(2)
          bc_coords(3) = 1.0D0 - bc_coords(1) - bc_coords(2)         

          if(minval(bc_coords) .ge. -eps) then
             id(iv) = itri
             cycle vecloop
          end if
          
       end do triangleloop
    end do vecloop

    return
       
  end function search_tr_vec
     
     
end module search_module

