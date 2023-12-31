!Method 1: Explicit backward difference
function explicit_bw (x, dt, N, u_IC, u_BC, c) result(u_np1) 
	implicit none
	
	real, dimension(1:), intent(in) :: x, u_IC
	real, intent(in) 				:: dt, u_BC, c
	integer, intent(in) 			:: N
	real, dimension(size(x))		:: u_n,u_np1
	
	integer :: idx, ndx !spacial, time indicies
	real :: dx !assumed uniform
	
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

!Method 2 Explicit Forward Difference. Unstable at this c
function explicit_fw (x, dt, N, u_IC, u_BC, c) result(u_np1) 
	implicit none
	
	real, dimension(1:), intent(in) :: x, u_IC
	real, intent(in) 				:: dt, u_BC, c
	integer, intent(in) 			:: N
	real, dimension(size(x))		:: u_n,u_np1
	
	integer :: idx, ndx !spacial, time indicies
	real :: dx !assumed uniform
	
	u_n = u_IC
	dx = (x(2)-x(1))
	do ndx = 1,N
		u_n(1) = u_BC 
		u_np1(1) = u_BC 
		do idx = 2,size(x)-1
			u_np1(idx) = u_n(idx)-(c*dt)/dx*(u_n(idx+1)-u_n(idx))
		end do
		!treat final node as bw difference
		idx = size(x) 
		u_np1(idx) = u_n(idx)-(c*dt)/dx*(u_n(idx)-u_n(idx-1))
		u_n = u_np1
	end do		
end function explicit_fw

!Method 3 Explicit Central Difference. *always* unstable
function explicit_cent (x, dt, N, u_IC, u_BC, c) result(u_np1) 
	implicit none
	
	real, dimension(1:), intent(in) :: x, u_IC
	real, intent(in) 				:: dt, u_BC, c
	integer, intent(in) 			:: N
	real, dimension(size(x))		:: u_n,u_np1
	
	integer :: idx, ndx !spacial, time indicies
	real :: dx !assumed uniform
	
	u_n = u_IC
	dx = (x(2)-x(1))
	do ndx = 1,N
		u_n(1) = u_BC 
		u_np1(1) = u_BC 
		do idx = 2,size(x)-1
			u_np1(idx) = u_n(idx)-(c*dt)/(2*dx)*(u_n(idx+1)-u_n(idx-1))
		end do
		!treat final node as bw difference
		idx = size(x) 
		u_np1(idx) = u_n(idx)-(c*dt)/dx*(u_n(idx)-u_n(idx-1))
		u_n = u_np1
	end do		
end function explicit_cent

!Method 6 Explicit Lax Scheme
function explicit_lax (x, dt, N, u_IC, u_BC, c) result(u_np1) 
	implicit none
	
	real, dimension(1:), intent(in) :: x, u_IC
	real, intent(in) 				:: dt, u_BC, c
	integer, intent(in) 			:: N
	real, dimension(size(x))		:: u_n,u_np1
	
	integer :: idx, ndx !spacial, time indicies
	real :: dx !assumed uniform
	
	u_n = u_IC
	dx = (x(2)-x(1))
	do ndx = 1,N
		u_n(1) = u_BC 
		u_np1(1) = u_BC 
		do idx = 2,size(x)-1
			u_np1(idx) = (u_n(idx-1)+u_n(idx+1))/2d0 -(c*dt)/(2*dx)*(u_n(idx+1)-u_n(idx-1))
		end do
		!treat final node as bw difference
		idx = size(x) 
		u_np1(idx) = u_n(idx)-(c*dt)/dx*(u_n(idx)-u_n(idx-1))
		u_n = u_np1
	end do		
end function explicit_lax

!Method 7 Explicit Lax-Wendroff 
function explicit_lax_wend(x, dt, N, u_IC, u_BC, c) result(u_np1) 
	implicit none
	
	real, dimension(1:), intent(in) :: x, u_IC
	real, intent(in) 				:: dt, u_BC, c
	integer, intent(in) 			:: N
	real, dimension(size(x))		:: u_n,u_np1
	
	integer :: idx, ndx !spacial, time indicies
	real :: dx !assumed uniform
	
	u_n = u_IC
	dx = (x(2)-x(1))
	do ndx = 1,N
		u_n(1) = u_BC 
		u_np1(1) = u_BC 
		do idx = 2,size(x)-1
			u_np1(idx) = u_n(idx)-(c*dt)/(2d0*dx)*(u_n(idx+1)-u_n(idx-1)) + &
			(0.5d0*((c**2d0)*(dt**2d0))/(dx**2d0))*(u_n(idx+1)-2d0*u_n(idx)+u_n(idx-1))
		end do
		!treat final node as bw difference
		idx = size(x) 
		u_np1(idx) = u_n(idx)-(c*dt)/dx*(u_n(idx)-u_n(idx-1))
		u_n = u_np1
	end do		
end function explicit_lax_wend

!Method 8 Explicit MacCormack
function explicit_macc(x, dt, N, u_IC, u_BC, c) result(u_np1) 
	implicit none
	
	real, dimension(1:), intent(in) :: x, u_IC
	real, intent(in) 				:: dt, u_BC, c
	integer, intent(in) 			:: N
	real, dimension(size(x))		:: u_n,u_np1,ub_np1
	
	integer :: idx, ndx !spacial, time indicies
	real :: dx !assumed uniform
	
	u_n = u_IC
	dx = (x(2)-x(1))
	do ndx = 1,N
		u_n(1) = u_BC 
		u_np1(1) = u_BC 
		ub_np1(1) = u_BC !need to do this? 
		do idx = 2,size(x)-1 
			ub_np1(idx) = u_n(idx)-c*dt/dx*(u_n(idx+1)-u_n(idx))
		end do
		do idx = 2,size(x)-1 
			u_np1(idx) = 0.5d0*(u_n(idx)+ub_np1(idx)-c*dt/dx*(ub_np1(idx)-ub_np1(idx-1)))
		end do
		!treat final node as bw difference
		idx = size(x) 
		u_np1(idx) = u_n(idx)-(c*dt)/dx*(u_n(idx)-u_n(idx-1))
		u_n = u_np1
	end do		
end function explicit_macc

!Method 9 Explicit Jameson 
function explicit_james(x, dt, N, u_IC, u_BC, c) result(u_np1) 
	implicit none
	
	real, dimension(1:), intent(in) :: x, u_IC
	real, intent(in) 				:: dt, u_BC, c
	integer, intent(in) 			:: N
	real, dimension(size(x))		:: u_n,u_np1, u_k, u_km1
	
	integer :: idx, ndx, kdx !spacial, time, step indicies
	real :: dx !assumed uniform
	real :: alp_k
	
	u_n = u_IC
	dx = (x(2)-x(1))
	
	do ndx = 1,N
		u_n(1) = u_BC
		u_np1(1) = u_BC
		u_k(1) = u_BC
		u_km1(1) = u_BC
		
		u_km1 = u_n !k=0
		do kdx = 1,4
			u_km1(1) = u_BC
			u_k(1) = u_BC
			alp_k = 1d0/(5d0-real(kdx))
			do idx = 2,size(x)-1
				u_k(idx) = u_n(idx) - alp_k*c*dt/(2d0*dx)*(u_km1(idx+1)-u_km1(idx-1))
			end do
			!treat final node as bw difference
			idx = size(x) 
			u_k(idx) = u_km1(idx)-(c*dt)/dx*(u_km1(idx)-u_km1(idx-1))
			u_km1 = u_k
		end do
		!treat final node as bw difference
		idx = size(x) 
		u_np1(idx) = u_n(idx)-(c*dt)/dx*(u_n(idx)-u_n(idx-1))
		u_np1 = u_k
		u_n = u_np1
	end do
	
end function explicit_james

!Method 10 Explicit Warming-Beam
function explicit_warm_beam(x, dt, N, u_IC, u_BC, c) result(u_np1) 
	implicit none
	
	real, dimension(1:), intent(in) :: x, u_IC
	real, intent(in) 				:: dt, u_BC, c
	integer, intent(in) 			:: N
	real, dimension(size(x))		:: u_n,u_np1,u_nph
	
	integer :: idx, ndx !spacial, time indicies
	real :: dx !assumed uniform
	
	u_n = u_IC
	dx = (x(2)-x(1))
	
	do ndx =  1,N
		u_n(1) = u_BC 
		u_np1(1) = u_BC
		u_nph(1) = u_BC
		do idx = 2,size(x)
			u_nph(idx) = u_n(idx) - c*dt/(2d0*dx)*(u_n(idx)-u_n(idx-1))
		end do
		do idx = 3, size(x) !may be able to do in one loop, but for simplicity
			u_np1(idx) = u_n(idx) - &
			c*dt/(dx)*((u_nph(idx)-u_nph(idx-1))+0.5d0*(u_n(idx)-2d0*u_n(idx-1)+u_n(idx-2)))
		end do
		
		idx = 2 !treat node 2 as bw difference (stencil too big, no ghost)
		u_np1(idx) = u_n(idx)-(c*dt)/dx*(u_n(idx)-u_n(idx-1))
		u_n = u_np1
	end do
end function explicit_warm_beam

!Method 11 Explicit Upwind
function explicit_upwind(x, dt, N, u_IC, u_BC, c) result(u_np1) 
	implicit none
	
	real, dimension(1:), intent(in) :: x, u_IC
	real, intent(in) 				:: dt, u_BC, c
	integer, intent(in) 			:: N
	real, dimension(size(x))		:: u_n,u_np1,fiph_n !prob dont need to store f
	
	integer :: idx, ndx !spacial, time indicies
	real :: dx !assumed uniform
	
	u_n = u_IC
	dx = (x(2)-x(1))
	
	
	do ndx = 1,N
		u_n(1) = u_BC 
		u_np1(1) = u_BC 
		do idx = 1,size(x)
			fiph_n(idx) = 0.5d0*c*(u_n(idx)+u_n(idx+1))-0.5d0*abs(c)*(u_n(idx+1)-u_n(idx))
		end do
		do idx = 2,size(x)
			u_np1(idx) = u_n(idx)-(dt/dx)*(fiph_n(idx)-fiph_n(idx-1))
		end do
		u_n = u_np1
	end do
	
		!treat final node as bw difference
		!idx = size(x) 
		!u_np1(idx) = u_n(idx)-(c*dt)/dx*(u_n(idx)-u_n(idx-1))		
end function explicit_upwind

subroutine output_result(path, x, u_IC, u)
	implicit none 
	character(*), intent(in) :: path
	real, dimension(1:), intent(in) :: x
	real, dimension(1:), intent(in)	:: u_IC
	real, dimension(1:), intent(in)	:: u
	
	!character, allocatable, intent(in) :: outpath(:)
	integer :: io
	integer :: idx
	integer :: I
	
	I = size(x)
	!allocate(outpath(size(path)))
	!outpath = path
	
	open (newunit=io, file=path, status="replace", action="write")
	write(io,*) I
	do idx=1, I
		write(io, *) x(idx),"  ", u_IC(idx),"  ", u(idx)
	end do
end subroutine output_result

program wave_eqn_ex
	implicit none
	
	!io
    integer :: io
    
	
	integer, parameter :: I = 41 !number mesh points
	integer, parameter :: N = 10 !number iterations
	real, parameter :: CFL = 0.9d0
	real, parameter :: c = 1.0d0
	real, parameter, dimension(2) :: xspan = [0d0, 2d0]
	real, parameter :: u0 = 1d0 !left BC
	
	real :: x(I)
	real :: dx
	real :: dt, dt2,dt3
	
	real,dimension(I) :: u
	real :: u_IC(I),u_exact1(I),u_exact2(I)
	
	interface !is this really the best way to do this? Lots of repeated code.
		function explicit_bw(x, dt, N, u_IC, u_BC, c) result(u_np1)
			real, dimension(1:), intent(in) :: x, u_IC
			real, intent(in) 				:: dt, u_BC, c
			integer, intent(in) 			:: N
			real, dimension(size(x))		:: u_np1
		end function
		function explicit_fw(x, dt, N, u_IC, u_BC, c) result(u_np1)
			real, dimension(1:), intent(in) :: x, u_IC
			real, intent(in) 				:: dt, u_BC, c
			integer, intent(in) 			:: N
			real, dimension(size(x))		:: u_np1
		end function
		function explicit_cent(x, dt, N, u_IC, u_BC, c) result(u_np1)
			real, dimension(1:), intent(in) :: x, u_IC
			real, intent(in) 				:: dt, u_BC, c
			integer, intent(in) 			:: N
			real, dimension(size(x))		:: u_np1
		end function
		function explicit_lax(x, dt, N, u_IC, u_BC, c) result(u_np1)
			real, dimension(1:), intent(in) :: x, u_IC
			real, intent(in) 				:: dt, u_BC, c
			integer, intent(in) 			:: N
			real, dimension(size(x))		:: u_np1
		end function
		function explicit_lax_wend(x, dt, N, u_IC, u_BC, c) result(u_np1)
			real, dimension(1:), intent(in) :: x, u_IC
			real, intent(in) 				:: dt, u_BC, c
			integer, intent(in) 			:: N
			real, dimension(size(x))		:: u_np1
		end function
		function explicit_macc(x, dt, N, u_IC, u_BC, c) result(u_np1)
			real, dimension(1:), intent(in) :: x, u_IC
			real, intent(in) 				:: dt, u_BC, c
			integer, intent(in) 			:: N
			real, dimension(size(x))		:: u_np1
		end function
		function explicit_james(x, dt, N, u_IC, u_BC, c) result(u_np1)
			real, dimension(1:), intent(in) :: x, u_IC
			real, intent(in) 				:: dt, u_BC, c
			integer, intent(in) 			:: N
			real, dimension(size(x))		:: u_np1
		end function
		function explicit_warm_beam(x, dt, N, u_IC, u_BC, c) result(u_np1)
			real, dimension(1:), intent(in) :: x, u_IC
			real, intent(in) 				:: dt, u_BC, c
			integer, intent(in) 			:: N
			real, dimension(size(x))		:: u_np1
		end function
		function explicit_upwind(x, dt, N, u_IC, u_BC, c) result(u_np1)
			real, dimension(1:), intent(in) :: x, u_IC
			real, intent(in) 				:: dt, u_BC, c
			integer, intent(in) 			:: N
			real, dimension(size(x))		:: u_np1
		end function
		
		subroutine output_result(path, x, u_IC, u)
			character(*), intent(in) :: path
			real, dimension(1:), intent(in) :: x
			real, dimension(1:), intent(in)	:: u_IC
			real, dimension(1:), intent(in)	:: u
		end subroutine
	end interface
	
	integer :: idx
	
	dx = (xspan(2)-xspan(1))/(real(I-1))
	do idx = 1,I
		x(idx) = (idx-1)*dx
	end do
	
	dt = (CFL*dx)/abs(c)
	dt2 = (2.0d0*dx)/abs(c)
	dt3 = (1.0d0*dx)/abs(c) !problem 5.2
	
	!IC prob should be defined w/ parameters
	do idx = 1,I 
		if (x(idx)<=0.5) then
			u_IC(idx) = 1
		else
			u_IC(idx) = 0.5
		end if
	end do
	
	!exact solutions
	
	do idx = 1,I 
		if (x(idx)<=(0.5d0+c*N*dt)) then
			u_exact1(idx) = 1
		else
			u_exact1(idx) = 0.5
		end if
	end do
	
	do idx = 1,I 
		if (x(idx)<=(0.5d0+c*N*dt2)) then
			u_exact2(idx) = 1
		else
			u_exact2(idx) = 0.5
		end if
	end do
	
	call output_result("exact_soln1.dat",x,u_IC,u_exact1)
	call output_result("exact_soln2.dat",x,u_IC,u_exact2)
	
	!method 1
	u = explicit_bw(x, dt, N, u_IC, u0, c)
	call output_result("explicit_bw.dat", x, u_IC, u)
	
	!method 2
	u = explicit_fw(x, dt, N, u_IC, u0, c)
	call output_result("explicit_fw.dat", x, u_IC, u)
	
	!method 3
	u = explicit_cent(x, dt, N, u_IC, u0, c)
	call output_result("explicit_cent.dat", x, u_IC, u)
	
	!method 6
	u = explicit_lax(x, dt, N, u_IC, u0, c)
	call output_result("explicit_lax.dat", x, u_IC, u)
	
	!method 7
	u = explicit_lax_wend(x, dt, N, u_IC, u0, c)
	call output_result("explicit_lax_wend.dat", x, u_IC, u)
	
	!method 8
	u = explicit_macc(x, dt, N, u_IC, u0, c)
	call output_result("explicit_macc.dat", x, u_IC, u)
	
	!method 9
	u = explicit_james(x, dt, N, u_IC, u0, c)
	call output_result("explicit_james.dat", x, u_IC, u)
	u = explicit_james(x, dt2, N, u_IC, u0, c) !TODO not quite right? (odd, 0.9 looks good)
	call output_result("explicit_james_CFL_2.0.dat", x, u_IC, u)
	
	!method 10
	u = explicit_warm_beam(x, dt, N, u_IC, u0, c)
	call output_result("explicit_warm_beam.dat", x, u_IC, u)
	u = explicit_warm_beam(x, dt2, N, u_IC, u0, c)
	call output_result("explicit_warm_beam_CFL_2.0.dat", x, u_IC, u)
	
	!method 11
	u = explicit_upwind(x, dt, N, u_IC, u0, c)
	call output_result("explicit_upwind.dat", x, u_IC, u)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Problem 5.2
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	!TODO: off by one? should be exact soln
	!method 1
	u = explicit_bw(x, dt3, N, u_IC, u0, c)
	call output_result("perfect_shift_bw.dat", x, u_IC, u)
	
	!method 6
	u = explicit_lax(x, dt3, N, u_IC, u0, c)
	call output_result("perfect_shift_lax.dat", x, u_IC, u)
	
	!method 7
	u = explicit_lax_wend(x, dt3, N, u_IC, u0, c)
	call output_result("perfect_shift_lax_wend.dat", x, u_IC, u)
	
	!method 8
	u = explicit_macc(x, dt3, N, u_IC, u0, c)
	call output_result("perfect_shift_macc.dat", x, u_IC, u)
	
	!method 10
	u = explicit_warm_beam(x, dt3, N, u_IC, u0, c)
	call output_result("perfect_shift_warm_beam.dat", x, u_IC, u)
	
	
end program wave_eqn_ex