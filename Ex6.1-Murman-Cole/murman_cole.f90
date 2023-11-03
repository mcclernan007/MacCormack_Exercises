! subroutine read_grid_file(filepath,xy)
	! implicit none
	! integer, parameter :: kint = SELECTED_INT_KIND(16) !general int, 16 digits
	! integer, parameter :: kflt = SELECTED_REAL_KIND(9,10) !general float, 9 digits, max val 1e10
	! integer, parameter :: kdub = SELECTED_REAL_KIND(20,200) !general "double", 20 digits,
	! integer :: io
	! logical :: exists
	! integer(kint) :: Npts
	! character(:), allocatable, intent(in) :: filepath
	! real(kdub), allocatable, intent(inout) :: xy(:,:)
	
	! integer(kint) :: idx
	! inquire(file=filepath, exist=exists)
	! if (exists) then
		! open (newunit=io, file=filepath, action="read")
		! read(io,*) Npts
		! allocate(xy(Npts,2))
		
		! do idx=1, Npts, 1
			! read(io,*) xy(idx,1), xy(idx,2)
		! end do
	  
		! close(io)
	! end if
! end subroutine read_grid_file

!interface 
!	real(kdub, allocatable(:,:) subroutine read_grid_file(filepath,xy):
!		real(kdub), allocatable, intent(inout) :: xy(:,:)
!		character(:), allocatable, intent(in) :: filepath

!end interface

program murman_cole
	implicit none
	 
    integer, parameter :: kint = SELECTED_INT_KIND(16) !general int, 16 digits
	integer, parameter :: kflt = SELECTED_REAL_KIND(9,10) !general float, 9 digits, max val 1e10
	integer, parameter :: kdub = SELECTED_REAL_KIND(20,200) !general "double", 20 digits,
    
!function return types
	real(kdub):: read_grid_file
	
	character(:), allocatable :: filepath
	real(kdub), allocatable :: xy(:,:)
	
	integer(kint) :: idx
	integer :: io
	logical :: exists
	integer(kint) :: Npts
	
	
	filepath = "./grid.txt"
	!call read_grid_file(grid_path,xy)
	
	
	inquire(file=filepath, exist=exists)
	if (exists) then
		open (newunit=io, file=filepath, action="read")
		read(io,*) Npts
		allocate(xy(Npts,2))
		
		do idx=1, Npts, 1
			read(io,*) xy(idx,1), xy(idx,2)
		end do
	  
		close(io)
	end if
	
	do idx = 1,size(xy,1),1
		print *, xy(idx,1), xy(idx,2)
	end do
	
end program murman_cole