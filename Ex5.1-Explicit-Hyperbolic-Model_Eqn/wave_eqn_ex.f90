function explicit_bw (x, dt, N, u_IC, u_BC, c) result(u_np1)
	implicit none
	
	real, dimension(1:), intent(in) :: x(:)
	real, intent(in) 				:: dt
	integer, intent(in) 			:: N
	real, dimension(1:), intent(in)	:: u_IC
	real, intent(in)				:: u_BC
	real, intent(in)				:: c
	
	
	integer :: idx !spacial index
	integer :: ndx !time index
	
	real :: dx !assumed uniform
	
	real, dimension(size(x)) :: u_n
	real, dimension(size(x)) :: u_np1
	
	u_n = u_IC
	dx = (x(2)-x(1))
	do ndx = 1,N
		u_n(1) = u_BC 
		u_np1(1) = u_BC 
		do idx = 2,size(x)
			u_np1(idx) = u_n(idx)-(c*dt)/dx*(u_n(idx)-u_n(idx-1))
		end do
		u_n = u_np1
	end do	
	
end function explicit_bw

program wave_eqn_ex
	implicit none
	
	!io
    integer :: io
    character(:), allocatable :: outPath
	
	integer, parameter :: I = 41 !number mesh points
	integer, parameter :: N = 10 !number iterations
	real, parameter :: CFL = 0.9d0
	real, parameter :: c = 1.0d0
	real, parameter, dimension(2) :: xspan = [0d0, 2d0]
	real, parameter :: u0 = 1d0 !left BC
	
	real :: x(I)
	real :: dx
	real :: dt
	
	real,dimension(I) :: u
	real :: u_IC(I) 
	
	interface
		function explicit_bw(x, dt, N, u_IC, u_BC, c) result(u_np1)
			real, dimension(1:), intent(in) :: x(:)
			real, intent(in) 				:: dt
			integer, intent(in) 			:: N
			real, dimension(1:), intent(in)	:: u_IC
			real, intent(in)				:: u_BC
			real, intent(in)				:: c
			real, dimension(size(x)) :: u_np1
		end function
	end interface
	
	integer :: idx
	
	dx = (xspan(2)-xspan(1))/(real(I-1))
	do idx = 1,I
		x(idx) = (idx-1)*dx
	end do
	
	dt = (CFL*dx)/abs(c)
	
	do idx = 1,I !IC prob should be defined w/ parameters
		if (x(idx)<=0.5) then
			u_IC(idx) = 1
		else
			u_IC(idx) = 0.5
		end if
	end do
	
	u = explicit_bw(x, dt, N, u_IC, u0, c)!method 1
	
	outPath = "explicit_bw.dat" !functionalize
	open (newunit=io, file=outPath, status="replace", action="write")
	write(io,*) I
	do idx=1, I
		write(io, *) x(idx),"  ", u_IC(idx),"  ", u(idx)
	end do
	 
	
	!u = explicit_fw(x, u_IC, dt, N)!method 2
	!u = explicit_cent(x, u_IC, dt, N)!method 3
	!u = explicit_lax(x, u_IC, dt, N)!method 6
	!u = explicit_lax_wend(x, u_IC, dt, N)!method 7
	!u = explicit_macc(x, u_IC, dt, N)!method 8
	!u = explicit_james(x, u_IC, dt, N)!method 9
	!u = explicit_warm_beam(x, u_IC, dt, N)!method 10
	!u = explicit_upwind(x, u_IC, dt, N)!method 11
	
end program wave_eqn_ex