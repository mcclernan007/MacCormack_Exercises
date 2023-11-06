program murman_cole
	implicit none
	 
	integer, parameter :: kint = SELECTED_INT_KIND(16) !general int, 16 digits
	integer, parameter :: kflt = SELECTED_REAL_KIND(9,10) !general float, 9 digits, max val 1e10
	integer, parameter :: kdub = SELECTED_REAL_KIND(20,200) !general "double", 20 digits,
	
	integer :: io
	logical :: exists
	
	character(:), allocatable :: grid_path
	real(kdub) :: M
	real(kdub) :: curX
	integer(kint) :: nstop
	
	real(kdub), allocatable :: nodes_xy(:,:) !dimensions=M*Nx2 format i1 j1; i1 j2;... i1 jn; i2 j1;... im jn
	real(kdub), allocatable :: grid_x(:,:) !dimensions M*N
	real(kdub), allocatable :: grid_y(:,:) !dimensions M*N
	real(kdub), allocatable :: A(:,:)      !dimensions M*N
	
	integer(kint) :: idx
	integer(kint) :: jdx
	integer(kint) :: kdx
	integer(kint) :: Npts !prob wasteful to save this too. 
	integer(kint) :: Nxpts
	integer(kint) :: Nypts
	
	
	M = 0.735_kdub
	!M = 0.908_kdub
	nstop = 400_kint
	grid_path = "./grid.txt"
	
	
	!should do in function or subroutine, but don't quite know how to return an allocatable array. Probably a syntatic problem
	!call read_grid_file(grid_path,nodes_xy)
	inquire(file=grid_path, exist=exists)
	if (exists) then
		open (newunit=io, file=grid_path, action="read")
		read(io,*) Npts, Nxpts, Nypts
		allocate(nodes_xy(Npts,2))
		allocate(grid_x(Nxpts,Nypts))
		allocate(grid_y(Nxpts,Nypts))
		kdx = 1
		do idx=1, Nxpts, 1
			do jdx=1, Nypts, 1
				read(io,*) grid_x(idx,jdx),grid_y(idx,jdx)
				nodes_xy(kdx,1) = grid_x(idx,jdx)
				nodes_xy(kdx,2) = grid_y(idx,jdx)
				kdx = kdx+1
			end do
		end do
	  
		close(io)
		else
			print *, "could not find a grid at path: " , grid_path
			error stop
	end if
	allocate(A(Nxpts,Nypts))
	allocate(mu(Nxpts,Nypts))
	allocate(phi_n(Nxpts,Nypts))
	allocate(phi_np1(Nxpts,Nypts))
	
end program murman_cole