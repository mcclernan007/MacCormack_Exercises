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

subroutine output_result(path, xy, phi_n,u,v,I,J)
	implicit none 
	character(*), intent(in) :: path
	integer, intent(in) :: I,J
	real, dimension(2,I,J), intent(in) :: xy
	real, dimension(I,J), intent(in) :: phi_n,u,v
	!character, allocatable, intent(in) :: outpath(:)
	integer :: io
	integer :: idx,jdx
	
	open (newunit=io, file=path, status="replace", action="write")
	write(io,*) I*J,I,J
	do idx=1, I
		do jdx=1,J
			write(io, *) xy(1,idx,jdx),"  ", xy(2,idx,jdx),"  ", phi_n(idx,jdx),"  ",u(idx,jdx),"  ",v(idx,jdx)
		end do
	end do
	close(io)
end subroutine output_result

subroutine output_cp(path,xy,u,v,I,J,M,gama,P_inf,rho_inf)
	implicit none 
	character(*), intent(in) 			:: path
	integer, intent(in) 				:: I,J
	real, dimension(2,I,J), intent(in)	:: xy
	real, dimension(I,J), intent(in)	:: u,v
	real, intent(in)					:: M,gama,P_inf,rho_inf
	
	integer 			:: io
	integer 			:: idx,jdx
	real				:: P
	real, dimension(I)	:: xCp, Cp
	
	!allocate(xCp(I))
	!allocate(Cp(I))
	
	do idx = 1,I
		xCp(idx) = xy(1,idx,1)
		P = P_inf*((1d0-(gama-1d0)/(2d0)*(M**2d0)*((u(idx,1)**2+v(idx,1)**2/(M*sqrt(gama*P_inf/rho_inf))-1d0)))**((gama-1d0)/gama))
		Cp(idx) = (P_inf-P)/(0.5d0*rho_inf* (M*sqrt(gama*P_inf/rho_inf))**2d0)
	end do
	
	open (newunit=io, file=path, status="replace", action="write")
	write(io,*) I
	do idx=1, I
		write(io, *) xCp(idx),"  ", Cp(idx)
	end do
	close(io)
end subroutine output_cp

function get_IC(xy,I,J,c,th,M,gama,P_inf,rho_inf) result(phi_IC)
	implicit none
	integer, intent(in)					:: I,J
	real, dimension(2,I,J), intent(in) 	:: xy
	real, intent(in)					:: c, th, M, gama,P_inf,rho_inf
	
	real, dimension(I,J)				:: phi_IC
	real, dimension(I,J)				:: x,y
	real								:: A, r, dphidy, xb
	integer								:: idx
	
	x = xy(1,:,:)
	y = xy(2,:,:)
	
	!assign initial	conditon
	phi_IC(:,:) = M*sqrt(gama*P_inf/rho_inf)*x(:,:)
	A = 1d0-M**2d0
	r = (c**2d0+th**2d0)/(4d0*th) !computed radius of circle for airfoil for convinienc
	
	!Set BC ICs
	!bottom is either dphi/dy=0 (symmetry) or dphi/dy = V_inf*(dy/dx)_body (slip no pen)
	do idx=1,I
		if (x(idx,1)<0d0 .or. x(idx,1)>c) then	
			dphidy = 0d0
		else
			xb = r*cos(atan2((r-0.5e0*th),(x(idx,1)-c/2e0)))+0.5d0*c!body x (x is chord line)
			dphidy = M*sqrt(gama*P_inf/rho_inf) * (xb-c/2d0)/(r*sqrt(1d0-((xb**2d0-c/2d0)/(r))**2))
		end if
		phi_IC(idx,1) = phi_IC(idx,2) - dphidy * (y(idx,2)-y(idx,1)) !bottom
	end do
end function get_IC

subroutine get_uv(xy, I,J, phi_n, u,v)
	implicit none
	
	integer, intent(in)					:: I,J
	real, dimension(2,I,J), intent(in) 	:: xy
	real, dimension(I,J), intent(in)	:: phi_n
	
	real, dimension(I,J), intent(out)	:: u,v
	
	real, dimension(I,J)				:: x,y
	integer								:: idx,jdx
	
	x = xy(1,:,:)
	y = xy(2,:,:)
	
	u(:,:) = 0d0
	v(:,:) = 0d0 
	
	!main volume
	do idx = 2,I-1 
		do jdx = 2,J-1
			u(idx,jdx) = (phi_n(idx+1,jdx)-phi_n(idx-1,jdx))/(x(idx+1,jdx)-x(idx-1,jdx))
			v(idx,jdx) = (phi_n(idx,jdx+1)-phi_n(idx,jdx-1))/(y(idx,jdx+1)-x(idx,jdx-1))
		end do
	end do
	
	!edges
	do idx = 2,I-1
		jdx = 1
		u(idx,jdx) = (phi_n(idx+1,jdx)-phi_n(idx-1,jdx))/(x(idx+1,jdx)-x(idx-1,jdx))
		jdx = J
		u(idx,jdx) = (phi_n(idx+1,jdx)-phi_n(idx-1,jdx))/(x(idx+1,jdx)-x(idx-1,jdx))
	end do
	
	do jdx = 2,J-1
		idx = 1
		v(idx,jdx) = (phi_n(idx,jdx+1)-phi_n(idx,jdx))/(y(idx,jdx+1)-y(idx,jdx))
		idx = I
		v(idx,jdx) = (phi_n(idx,jdx+1)-phi_n(idx,jdx-1))/(y(idx,jdx+1)-x(idx,jdx-1))
	end do
	
	!corners
	jdx = 1
	u(1,jdx) = (phi_n(2,jdx)-phi_n(1,jdx))/(x(2,jdx)-x(1,jdx))
	u(I,jdx) = (phi_n(I,jdx)-phi_n(I-1,jdx))/(x(I,jdx)-x(I-1,jdx))
	jdx = J
	u(1,jdx) = (phi_n(2,jdx)-phi_n(1,jdx))/(x(2,jdx)-x(1,jdx))
	u(I,jdx) = (phi_n(I,jdx)-phi_n(I-1,jdx))/(x(I,jdx)-x(I-1,jdx))
		
	idx = 1
	v(idx,1) = (phi_n(idx,2)-phi_n(idx,1))/(y(idx,2)-y(idx,1))
	v(idx,J) = (phi_n(idx,J)-phi_n(idx,J-1))/(y(idx,J)-y(idx,J-1))
	idx = I
	v(idx,1) = (phi_n(idx,2)-phi_n(idx,1))/(y(idx,2)-y(idx,1))
	v(idx,J) = (phi_n(idx,J)-phi_n(idx,J-1))/(y(idx,J)-y(idx,J-1))

end subroutine get_uv

function point_jacobi(xy,phi_IC, nstop, I,J, c,th,M,gama,P_inf,rho_inf) result (phi_n)
	implicit none
	
	integer, intent(in)					:: nstop,I,J
	real, dimension(2,I,J), intent(in) 	:: xy
	real, dimension(I,J), intent(in)	:: phi_IC
	real, intent(in)					:: c, th, M, gama,P_inf,rho_inf
	
	real, allocatable 	:: phi_n(:,:), phi_np1(:,:),x(:,:),y(:,:)
	integer				:: ndx,idx,jdx,io
	real				:: A, Bi, dtxi, dtyi, dphidy,  r, xb 
	real				:: res, lres
	
	allocate(phi_n(I,J))
	allocate(phi_np1(I,J))
	allocate(x(I,J))
	allocate(y(I,J))
	x = xy(1,:,:)
	y = xy(2,:,:)
	
	A = 1d0-M**2d0
	r = (c**2d0+th**2d0)/(4d0*th)
	
	phi_n(:,:) = phi_IC(:,:)
	
	open (newunit=io, file="1-point_jacobi-resid.dat", status="replace", action="write")
	
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
			phi_n(1,jdx)   = M*sqrt(gama*P_inf/rho_inf)*x(1,jdx) !left
			phi_np1(1,jdx) = M*sqrt(gama*P_inf/rho_inf)*x(1,jdx) !left
			phi_n(I,jdx)   = M*sqrt(gama*P_inf/rho_inf)*x(I,jdx) !right
			phi_np1(I,jdx) = M*sqrt(gama*P_inf/rho_inf)*x(I,jdx) !right
		end do 
		
		do idx=1,I
			phi_n(idx,J)   = M*sqrt(gama*P_inf/rho_inf)*x(idx,J) !top 
			phi_np1(idx,J) = M*sqrt(gama*P_inf/rho_inf)*x(idx,J) !top 
		end do 
		
		!bottom is either dphi/dy=0 (symmetry) or dphi/dy = V_inf*(dy/dx)_body (slip no pen)
		do idx=1,I
			if (x(idx,1)<0d0 .or. x(idx,1)>c) then	
				dphidy = 0d0
			else
				xb = r*cos(atan2((r-0.5e0*th),(x(idx,1)-c/2e0)))+0.5d0*c!body x (x is chord line)
				dphidy = M*sqrt(gama*P_inf/rho_inf) * (xb-c/2d0)/(r*sqrt(1d0-((xb**2d0-c/2d0)/(r))**2))
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
end function point_jacobi

function point_gauss_seidel(xy,phi_IC, nstop, I,J, c,th,M,gama,P_inf,rho_inf) result (phi_n)
	implicit none
	
	integer, intent(in)					:: nstop,I,J
	real, dimension(2,I,J), intent(in) 	:: xy
	real, dimension(I,J), intent(in)	:: phi_IC
	real, intent(in)					:: c, th, M, gama,P_inf,rho_inf
	
	real, allocatable 	:: phi_n(:,:),x(:,:),y(:,:)
	integer				:: ndx,idx,jdx,io
	real				:: A, Bi, dtxi, dtyi, dphidy,  r, xb 
	real				:: res, lres
	
	allocate(phi_n(I,J))
	allocate(x(I,J))
	allocate(y(I,J))
	x = xy(1,:,:)
	y = xy(2,:,:)
	
	A = 1d0-M**2d0
	r = (c**2d0+th**2d0)/(4d0*th)
	
	phi_n(:,:) = phi_IC(:,:)
	
	open (newunit=io, file="2-point_gauss_sidel-resid.dat", status="replace", action="write")
	
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
				phi_n(idx,jdx) = 1d0/Bi*(A/dtxi*(phi_n(idx+1,jdx)/(x(idx+1,jdx)-x(idx,jdx))) + &
								  A/dtxi*(phi_n(idx-1,jdx)/(x(idx,jdx)-x(idx-1,jdx))) + &
								1d0/dtyi*(phi_n(idx,jdx+1)/(y(idx,jdx+1)-y(idx,jdx))) + &
								1d0/dtyi*(phi_n(idx,jdx-1)/(y(idx,jdx)-y(idx,jdx-1))))
			end do 
		end do
		
		!Set BCs
		!inlet,outlet, and top are const phi
		do jdx=1,J
			phi_n(1,jdx)   = M*sqrt(gama*P_inf/rho_inf)*x(1,jdx) !left
			phi_n(I,jdx)   = M*sqrt(gama*P_inf/rho_inf)*x(I,jdx) !right
		end do 
		
		do idx=1,I
			phi_n(idx,J)   = M*sqrt(gama*P_inf/rho_inf)*x(idx,J) !top  
		end do 
		
		!bottom is either dphi/dy=0 (symmetry) or dphi/dy = V_inf*(dy/dx)_body (slip no pen)
		do idx=1,I
			if (x(idx,1)<0d0 .or. x(idx,1)>c) then	
				dphidy = 0d0
			else
				xb = r*cos(atan2((r-0.5e0*th),(x(idx,1)-c/2e0)))+0.5d0*c!body x (x is chord line)
				dphidy = M*sqrt(gama*P_inf/rho_inf) * (xb-c/2d0)/(r*sqrt(1d0-((xb**2d0-c/2d0)/(r))**2))
			end if
			phi_n(idx,1)   = phi_n(idx,2)-dphidy*(x(idx,2)-x(idx,1)) !bottom
		end do 
		
		!compute and print residual, print, continue
		!res = norm2(phi_n-phi_np1)		
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
end function point_gauss_seidel

program elliptic
	implicit none
	
	real, allocatable 	:: phi_IC(:,:), phi_n(:,:),u(:,:),v(:,:)
	real, allocatable 	:: xCp(:),Cp(:)
	real, allocatable 	:: xy(:,:,:), uv(:,:,:) 
	
	integer				:: idx,jdx, I, J, ndx, io
	real				:: A, Bi, dtxi, dtyi, r, dphidy, xb, P! i is indicating index, not inv
	real				:: lres,res 
	
	real, parameter 	:: M = 0.5d0
	real, parameter 	:: P_inf = 1d0
	real, parameter 	:: rho_inf = 1d0
	real, parameter 	:: gama = 1.4d0
	real, parameter		:: c = 1.0d0 !may be able to interpret from exteranl grid. For simplicity...
	real, parameter		:: th = 0.06d0
	integer, parameter 	:: nstop = 1000
	
	interface 
		function read_grid_file(grid_path) result(grid_xy)
			character(*), intent(in)	:: grid_path
			real, allocatable 			:: grid_xy(:,:,:)
		end function read_grid_file
		function point_jacobi(xy,phi_IC, nstop, I,J, c,th,M,gama,P_inf,rho_inf) result(phi_n)
			integer, intent(in)					:: nstop,I,J
			real, dimension(2,I,J), intent(in) 	:: xy
			real, dimension(I,J), intent(in)	:: phi_IC
			real, intent(in)					:: c, th, M, gama,P_inf,rho_inf
			real, allocatable 	:: phi_n(:,:)
		end function point_jacobi
		function point_gauss_seidel(xy,phi_IC, nstop, I,J, c,th,M,gama,P_inf,rho_inf) result(phi_n)
			integer, intent(in)					:: nstop,I,J
			real, dimension(2,I,J), intent(in) 	:: xy
			real, dimension(I,J), intent(in)	:: phi_IC
			real, intent(in)					:: c, th, M, gama,P_inf,rho_inf
			real, allocatable 	:: phi_n(:,:)
		end function point_gauss_seidel
		function get_IC(xy,I,J,c,th,M,gama,P_inf,rho_inf) result(phi_IC)
			integer, intent(in)					:: I,J
			real, dimension(2,I,J), intent(in) 	:: xy
			real, intent(in)					:: c, th, M, gama,P_inf,rho_inf
			real, dimension(I,J)				:: phi_IC
		end function get_IC
	end interface
	
	xy = read_grid_file("grid.txt") 
	I = size(xy,2)
	J = size(xy,3)
	
	allocate(phi_IC(I,J))
	allocate(phi_n(I,J))
	allocate(u(I,J))
	allocate(v(I,J))
	
	
	phi_IC = get_IC(xy,I,J,c,th,M,gama,P_inf,rho_inf)
	call get_uv(xy, I,J, phi_IC, u,v)	
	call output_result("0-inital_conditions.dat", xy, phi_IC,u,v,I,J)
	
	!Method 1: point jacobi	
	phi_n = point_jacobi(xy,phi_IC, nstop, I,J, c,th,M,gama,P_inf,rho_inf)
	call get_uv(xy, I,J, phi_n, u,v)	
	call output_result("1-point_jacobi.dat", xy, phi_n,u,v,I,J)
	call output_cp("1-cp_point_jacobi.dat",xy,u,v,I,J,M,gama,P_inf,rho_inf)
	
	!Method 2: Point Gauss-Seidel
	phi_n = point_gauss_seidel(xy,phi_IC, nstop, I,J, c,th,M,gama,P_inf,rho_inf)
	call get_uv(xy, I,J, phi_n, u,v)	
	call output_result("2-point_gauss_siedel.dat", xy, phi_n,u,v,I,J)
	call output_cp("2-cp_point_gauss_siedel.dat",xy,u,v,I,J,M,gama,P_inf,rho_inf)
	!
	!
	
	
	
	
	
	
end program elliptic