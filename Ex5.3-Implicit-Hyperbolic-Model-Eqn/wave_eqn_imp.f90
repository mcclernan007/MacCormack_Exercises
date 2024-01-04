!Method 4,5: Implicit alpha method. 
!  alpha=1: Implicit Central
!  alpha=1/2: Crank-Nicholson
function implicit_alph (x, dt, N, u_IC, u_BC, c, alph) result(u_out) 
	implicit none

	real, dimension(1:), intent(in) :: x, u_IC
	real, intent(in) 				:: dt, u_BC, c, alph
	integer, intent(in) 			:: N
	real, dimension(size(x))		:: u_out,v!,f
	real, dimension(size(x)+1)		:: gam,u_np1,u_n
	!real, dimension(size(x),size(x)):: A
	
	real 							:: ai, bi, ci, ali, fi
	
	integer :: idx, jdx, ndx !spacial, time indicies
	integer :: I
	real :: dx !assumed uniform
	
	I = size(x)
	
	u_n = u_IC
	dx = (x(2)-x(1))
	
	do ndx = 1,N
		!step 1 (not stated, but ghost nodes?)
		u_np1(I+1) = 0d0
		v(I+1) = 0d0
		gam(I+1) = 0d0
		
		!step 2. tridiagonal inv
		idx = I
		ai = 1d0 + c*dt/dx
		bi = 0d0
		ci = -c*dt/dx
		ali = ai
		!calc, store
		
		gam(idx) = ci/ali
		v(idx) = ((u_n(idx))-bi*u_np1(idx+1))/ali
		do idx = I-1,2,-1
			!calc, dont store
			ai = 1d0
			bi = alph*(c*dt/(2d0*dx))
			ci = -bi
			ali = ai-bi*gam(idx+1)
			
			!calc, store
			gam(idx) = ci/ali
			fi = u_n(idx)-(1-alph)*c*dt/(2*dx)*(u_n(idx+1)-u_n(idx-1))
			v(idx) = (fi-bi*v(idx+1))/ali
		end do
		
		!calc, dont store
		idx = 1
		ai = 1
		bi = 0
		ci = 0
		ali = ai-bi*gam(idx+1)
		!calc, store
		gam(idx) = ci/ali
		v(idx) = (u_BC - ci*u_BC - bi*v(idx+1))/ali
		
		!step3 backsolve
		u_np1(1) = v(1)
		do idx=2,I
			u_np1(idx) = v(idx)-gam(idx)*u_np1(idx-1)
		end do
		u_n = u_np1
	end do
	u_out = u_n(1:I)
	
	!form A and f
	
	!do idx = 1,I
	!	do jdx = 1,I
	!		A(idx,jdx) = 0d0
	!	end do
	!end do
	
	!A(1,1) = 1
	!A(1,2) = -alph*(c*dt/(2*dx))
	!A(I,I) = 1
	!A(I,I-1) = alph*(c*dt/(2*dx))
	!do idx = 2,I-1
	!	A(idx,idx) = 1d0 !ai
	!	A(idx,idx-1) = alph*(c*dt/(2*dx))
	!	A(idx,idx+1) = -alph*(c*dt/(2*dx))
	!end do
	!
	!!Note f is RHS of Au = f. Not quite equivelent to f_i defined in maccormick
	!!ordering is odd in maccormick. [I to 1]?
	!f(1) = u_n(idx)-(1-alph)*c*dt/(2*dx)*(u_n(idx+1)-u_n(idx-1))+alph*(c*dt/(2*dx))*u_BC
	!f(I) = u_n(idx)
	!do idx = 2,I
	!	f(idx) = u_n(idx) - (1-alph)*c*dt/(2*dx)*(u_n(idx+1)-u_n(idx-1))
	!end do

	!from here would solve. not sure this explicit writing out of matricies is necc

	!u_np1(idx) = u_n(idx) - (1-alph)*c*dt/2/dx*(u_n(idx+1)-u_n(idx-1)) &
	! - alph*c*dt/2/dx*(u_np1(idx+1)-u_np1(idx-1))
	
	!Form matrix A !may not need to do explicitly (but dont have book infront of me)
	!good check of the thomas alg. Also would need to import lappack, which s a good exercise.
	
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

program wave_eqn_imp
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
	real :: dx, dt1, dt2
	
	interface
		function implicit_alph(x, dt, N, u_IC, u_BC, c, alph) result(u_n)
			real, dimension(1:), intent(in) :: x, u_IC
			real, intent(in) 				:: dt, u_BC, c, alph
			integer, intent(in) 			:: N
			real, dimension(size(x))		:: u_n
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
	
	!method 4 fully implicit, central (alph=1.0d0)
	u = implicit_alph(x, dt1, N, u_IC, u0, c, 1.0e0)
	call output_result("implicit_central.dat", x, u_IC, u)
	u = implicit_alph(x, dt2, N, u_IC, u0, c, 1.0e0)
	call output_result("implicit_central_CFL2.dat", x, u_IC, u)
	
	!method 5 crank-nicholson (implicit) (alph=0.5d0)
	u = implicit_alph(x, dt1, N, u_IC, u0, c, 0.5e0)
	call output_result("crank_nich.dat", x, u_IC, u)
	u = implicit_alph(x, dt2, N, u_IC, u0, c, 0.5e0)
	call output_result("crank_nich_CFL2.dat", x, u_IC, u)
	
end program wave_eqn_imp