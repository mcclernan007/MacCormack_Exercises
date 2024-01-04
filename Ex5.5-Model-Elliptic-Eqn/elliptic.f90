function read_grid_file(grid_path) result(grid_xy)
	implicit none
	character(*), intent(in)	:: grid_path
	real, allocatable		 	:: grid_xy(:,:,:)
	integer						:: Npts, Nxpts, Nypts, idx, jdx, io
	logical						:: exists
	
	inquire(file=grid_path, exist=exists)
	if (exists) then
		open (newunit=io, file=grid_path, action="read")
		read(io,*) Npts, Nxpts, Nypts
		allocate(grid_xy(2,Nxpts,Nypts))
		do idx=1, Nxpts, 1
			do jdx=1, Nypts, 1
				!read(io,*) grid_x(idx,jdx),grid_y(idx,jdx)
				read(io,*) grid_xy(1,idx,jdx), grid_xy(2,idx,jdx)
			end do
		end do
	  
		close(io)
		else
			print *, "could not find a grid at path: " , grid_path
			error stop
	end if
	
end function read_grid_file


program elliptic
	implicit none
	
	real, allocatable 	:: xy(:,:,:) !wierd to use 3 dimensional, but fuctional 
	real, allocatable 	:: phi_n(:,:), phi_np1(:,:),x(:,:),y(:,:)
	
	integer				:: idx,jdx, I, J, ndx
	real				:: A, Bi, dtxi, dtyi ! i is indicating index, not inv
	
	real, parameter 	:: M = 0.5d0
	real, parameter 	:: P =1d0
	real, parameter 	:: rho = 1d0
	real, parameter 	:: gama = 1.4d0
	integer, parameter 	:: nstop = 2
	
	interface 
		function read_grid_file(grid_path) result(grid_xy)
			character(*), intent(in)	:: grid_path
			real, allocatable 			:: grid_xy(:,:,:)
		end function read_grid_file
	end interface
	
	xy = read_grid_file("grid.txt") 
	I = size(xy,2)
	J = size(xy,3)
	
	allocate(phi_n(I,J))
	allocate(phi_np1(I,J))
	allocate(x(I,J))
	allocate(y(I,J))
	x = xy(1,:,:)
	y = xy(2,:,:)
	!deallocate(xy) !not optimal, but for indicial consistency
	
	!assign initial	conditon
	phi_n(:,:) = M*sqrt(gama*P/rho)*xy(1,:,:)
	
	A = 1-M**2
	!Method 1: point jacobi
	do ndx = 1,nstop
		do idx = 2,I-1
			do jdx = 2,J-1
				!nonuniform spacing makes this super messy
				dtxi = (x(idx+1,jdx) - x(idx-1,jdx))/2d0
				dtyi = (y(idx,jdx+1) - y(idx,jdx-1))/2d0
				Bi = A/dtxi*(1d0/(x(idx+1,jdx)-x(idx,jdx)) + 1d0/(x(idx,jdx)-x(idx-1,jdx))) + &
				   1d0/dtyi*(1d0/(y(idx,jdx+1)-y(idx,jdx)) + 1d0/(y(idx,jdx)-y(idx,jdx-1)))
				
				!B of index i, not inverse
				phi_np1 = 1d0/Bi*(A/dtxi*(phi_n(idx+1,jdx)/(x(idx+1,jdx)-x(idx,jdx))) + &
								  A/dtxi*(phi_n(idx-1,jdx)/(x(idx,jdx)-x(idx-1,jdx))) + &
								1d0/dtyi*(phi_n(idx,jdx+1)/(y(idx,jdx+1)-y(idx,jdx))) + &
								1d0/dtyi*(phi_n(idx,jdx-1)/(y(idx,jdx)-y(idx,jdx-1))))
			end do 
		end do
	end do
	
end program elliptic