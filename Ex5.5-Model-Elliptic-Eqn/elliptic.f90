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

subroutine output_result(path, x, y, phi_n,u,v,I,J)
	implicit none 
	character(*), intent(in) :: path
	integer, intent(in) :: I,J
	real, dimension(I,J), intent(in) :: x,y,phi_n,u,v
	!character, allocatable, intent(in) :: outpath(:)
	integer :: io
	integer :: idx,jdx
	
	open (newunit=io, file=path, status="replace", action="write")
	write(io,*) I*J,I,J
	do idx=1, I
		do jdx=1,J
			write(io, *) x(idx,jdx),"  ", y(idx,jdx),"  ", phi_n(idx,jdx),"  ",u(idx,jdx),"  ",v(idx,jdx)
		end do
	end do
	close(io)
end subroutine output_result

subroutine output_cp(path, x, Cp, I)
	implicit none 
	character(*), intent(in) :: path
	integer, intent(in) :: I
	real, dimension(I), intent(in) :: x,Cp
	!character, allocatable, intent(in) :: outpath(:)
	integer :: io
	integer :: idx,jdx
	
	open (newunit=io, file=path, status="replace", action="write")
	write(io,*) I
	do idx=1, I
		write(io, *) x(idx),"  ", Cp(idx)
	end do
	close(io)
end subroutine output_cp


program elliptic
	implicit none
	
	real, allocatable 	:: phi_n(:,:), phi_np1(:,:),x(:,:),y(:,:),u(:,:),v(:,:)
	real, allocatable 	:: xCp(:),Cp(:)
	real, allocatable 	:: xy(:,:,:) !wierd to use 3 dimensional, but fuctional 
	
	integer				:: idx,jdx, I, J, ndx, io
	real				:: A, Bi, dtxi, dtyi, r, dphidy, xb, P! i is indicating index, not inv
	real				:: lres,res 
	
	real, parameter 	:: M = 0.5d0
	real, parameter 	:: P_inf = 1d0
	real, parameter 	:: rho = 1d0
	real, parameter 	:: gama = 1.4d0
	real, parameter		:: c = 1.0d0 !may be able to interpret from exteranl grid. For simplicity...
	real, parameter		:: th = 0.06d0
	integer, parameter 	:: nstop = 1000
	
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
	
	!assign initial	conditon
	phi_n(:,:) = M*sqrt(gama*P_inf/rho)*x(:,:)
	A = 1d0-M**2d0
	r = (c**2d0+th**2d0)/(4d0*th) !computed radius of circle for airfoil for convinience
	
	res = 1d-30
	lres = 1d-30
	!Method 1: point jacobi	
	
	
	!Set BC ICs
	!bottom is either dphi/dy=0 (symmetry) or dphi/dy = V_inf*(dy/dx)_body (slip no pen)
	do idx=1,I
		if (x(idx,1)<0d0 .or. x(idx,1)>c) then	
			dphidy = 0d0
		else
			xb = r*cos(atan2((r-0.5e0*th),(x(idx,1)-c/2e0)))+0.5d0*c!body x (x is chord line)
			dphidy = M*sqrt(gama*P_inf/rho) * (xb-c/2d0)/(r*sqrt(1d0-((xb**2d0-c/2d0)/(r))**2))
		end if
		phi_n(idx,1) = phi_n(idx,2) - dphidy * (y(idx,2)-y(idx,1)) !bottom
	end do 
	
	allocate(u(I,J))
	allocate(v(I,J))
	u(:,:) = 0d0
	v(:,:) = 0d0 
	do idx = 2,I-1 !neglect bounds for now
		do jdx = 2,J-1
			u(idx,jdx) = (phi_n(idx+1,jdx)-phi_n(idx-1,jdx))/(x(idx+1,jdx)-x(idx-1,jdx))
			v(idx,jdx) = (phi_n(idx,jdx+1)-phi_n(idx,jdx-1))/(y(idx,jdx+1)-x(idx,jdx-1))
		end do
	end do
	
	call output_result("inital_conditions.dat", x, y, phi_n,u,v,I,J)
	
	open (newunit=io, file="resid.dat", status="replace", action="write")
	
	
	do ndx = 1,nstop
		
		!Main loops
		do idx = 2,I-1
			do jdx = 2,J-1
				!nonuniform spacing makes this super messy
				dtxi = (x(idx+1,jdx) - x(idx-1,jdx))/2d0
				dtyi = (y(idx,jdx+1) - y(idx,jdx-1))/2d0
				Bi = A/dtxi*(1d0/(x(idx+1,jdx)-x(idx,jdx)) + 1d0/(x(idx,jdx)-x(idx-1,jdx))) + &
				   1d0/dtyi*(1d0/(y(idx,jdx+1)-y(idx,jdx)) + 1d0/(y(idx,jdx)-y(idx,jdx-1)))
				
				!B of index i, not inverse
				phi_np1(idx,jdx) = 1d0/Bi*(A/dtxi*(phi_n(idx+1,jdx)/(x(idx+1,jdx)-x(idx,jdx))) + &
								  A/dtxi*(phi_n(idx-1,jdx)/(x(idx,jdx)-x(idx-1,jdx))) + &
								1d0/dtyi*(phi_n(idx,jdx+1)/(y(idx,jdx+1)-y(idx,jdx))) + &
								1d0/dtyi*(phi_n(idx,jdx-1)/(y(idx,jdx)-y(idx,jdx-1))))
			end do 
		end do
		
		!Set BCs
		!inlet,outlet, and top are const phi
		do jdx=1,J
			phi_n(1,jdx)   = M*sqrt(gama*P_inf/rho)*x(1,jdx) !left
			phi_np1(1,jdx) = M*sqrt(gama*P_inf/rho)*x(1,jdx) !left
			phi_n(I,jdx)   = M*sqrt(gama*P_inf/rho)*x(I,jdx) !right
			phi_np1(I,jdx) = M*sqrt(gama*P_inf/rho)*x(I,jdx) !right
		end do 
		
		do idx=1,I
			phi_n(idx,J)   = M*sqrt(gama*P_inf/rho)*x(idx,J) !top 
			phi_np1(idx,J) = M*sqrt(gama*P_inf/rho)*x(idx,J) !top 
		end do 
		
		!bottom is either dphi/dy=0 (symmetry) or dphi/dy = V_inf*(dy/dx)_body (slip no pen)
		do idx=1,I
			if (x(idx,1)<0d0 .or. x(idx,1)>c) then	
				dphidy = 0d0
			else
				xb = r*cos(atan2((r-0.5e0*th),(x(idx,1)-c/2e0)))+0.5d0*c!body x (x is chord line)
				dphidy = M*sqrt(gama*P_inf/rho) * (xb-c/2d0)/(r*sqrt(1d0-((xb**2d0-c/2d0)/(r))**2))
			end if
			phi_n(idx,1)   = phi_n(idx,2)-dphidy*(x(idx,2)-x(idx,1)) !bottom
			phi_np1(idx,1) = phi_n(idx,2)-dphidy*(x(idx,2)-x(idx,1)) !bottom
		end do 
		
		!compute and print residual, print, continue
		!res = norm2(phi_n-phi_np1)
		phi_n = phi_np1
		
		res = 1d-30
		do idx = 2,I-1
			do jdx = 2,J-1
				lres = abs(A*( &
				(((phi_n(idx+1,jdx) - phi_n(idx,jdx))/(x(idx+1,jdx)-x(idx,jdx))) - &
				 ((phi_n(idx,jdx) - phi_n(idx-1,jdx))/(x(idx,jdx) - x(idx-1,jdx)))) / &
				 (0.5d0*(x(idx+1,jdx)-x(idx-1,jdx)))) +  &
				 ((((phi_n(idx,jdx+1) - phi_n(idx,jdx))/(y(idx,jdx+1)-y(idx,jdx))) - &
				  ((phi_n(idx,jdx) - phi_n(idx,jdx-1))/(y(idx,jdx) - y(idx,jdx-1)))) / &
				 (0.5d0*(y(idx,jdx+1)-y(idx,jdx-1))))) !oof
				if (lres>res) then 
					!print*, "location of max",idx,jdx
					res=lres
				end if
			end do
		end do
		print *,ndx, res
		write(io,*) ndx, res
	end do
	close(io)
	
	!P_infostprocess, get u,v
	
	!allocate(u(I,J))
	!allocate(v(I,J))
	u(:,:) = 0d0
	v(:,:) = 0d0 
	do idx = 2,I-1 
		do jdx = 2,J-1
			u(idx,jdx) = (phi_n(idx+1,jdx)-phi_n(idx-1,jdx))/(x(idx+1,jdx)-x(idx-1,jdx))
			v(idx,jdx) = (phi_n(idx,jdx+1)-phi_n(idx,jdx-1))/(y(idx,jdx+1)-x(idx,jdx-1))
		end do
	end do
	
	do idx = 1,I
		jdx = 1
		u(idx,jdx) = (phi_n(idx+1,jdx)-phi_n(idx,jdx))/(x(idx+1,jdx)-x(idx,jdx))
		jdx = J
		u(idx,jdx) = (phi_n(idx,jdx)-phi_n(idx-1,jdx))/(x(idx,jdx)-x(idx-1,jdx))
	end do
	
	do jdx = 1,J
		idx = 1
		v(idx,jdx) = (phi_n(idx,jdx+1)-phi_n(idx,jdx))/(y(idx,jdx+1)-y(idx,jdx))
		idx = I
		v(idx,J) = (phi_n(idx,jdx)-phi_n(idx,jdx-1))/(y(idx,jdx)-y(idx,jdx-1))
	end do
	
	call output_result("point_jacobi.dat", x, y, phi_n,u,v,I,J)
	
	allocate(xCp(I))
	allocate(Cp(I))
	do idx = 1,I
		xCp(idx) = x(idx,1)
		P = P_inf*((1d0-(gama-1d0)/(2d0)*(M**2d0)*((u(idx,1)**2+v(idx,1)**2/(M*sqrt(gama*P_inf/rho))-1d0)))**((gama-1d0)/gama))
		Cp(idx) = (P_inf-P)/(0.5d0*rho* (M*sqrt(gama*P_inf/rho))**2d0)
	end do
	
	call output_cp("cp_point_jacobi.dat", xCp, Cp, I)
	
	
	
	
	
end program elliptic