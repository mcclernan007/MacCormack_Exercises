!Method 4,5: Implicit alpha method. 
!  alpha=1: Implicit Central
!  alpha=1/2: Crank-Nicholson
function implicit_alph (x, dt, N, u_IC, u_BC, c, alph) result(u_np1) 
	implicit none
	
	real, dimension(1:), intent(in) :: x, u_IC
	real, intent(in) 				:: dt, u_BC, c, alph
	integer, intent(in) 			:: N
	real, dimension(size(x))		:: u_n,u_np1
	
	integer :: idx, ndx !spacial, time indicies
	real :: dx !assumed uniform
	
	u_n = u_IC
	dx = (x(2)-x(1))
	
end function implicit_alph


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
	real, parameter :: CFL1 = 0.9d0
	real, parameter :: CFL2 = 2.0d0
	
	real, parameter :: c = 1.0d0
	real, parameter, dimension(2) :: xspan = [0d0, 2d0]
	real, parameter :: u0 = 1d0 !left BC
	
	real, dimension(I) :: x, u, u_IC, u_exact1, u_exact2
	real :: dx, dt1, dt2,
	
	interface !is this really the best way to do this? Lots of repeated code.
		function implicit_alph(x, dt, N, u_IC, u_BC, c, alph) result(u_np1)
			real, dimension(1:), intent(in) :: x, u_IC
			real, intent(in) 				:: dt, u_BC, c, alph
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
	
	dt1 = (CFL1*dx)/abs(c)
	dt2 = (CFL2*dx)/abs(c)
	
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
		if (x(idx)<=(0.5d0+c*N*dt1)) then
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
	
	!method 4 fully implicit, central
	u = explicit_alpha(x, dt, N, u_IC, u0, c)
	call output_result("implicit_central.dat", x, u_IC, u, 1.0d0)
	
	!method 5 crank-nicholson (implicit)
	u = explicit_fw(x, dt, N, u_IC, u0, c)
	call output_result("crank_nich.dat", x, u_IC, u, 0.5d0)
	
end program wave_eqn_ex