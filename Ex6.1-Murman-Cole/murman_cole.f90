function read_grid_file(grid_path) result(grid_xy)
    implicit none
    character(*), intent(in)    :: grid_path
    real(8), allocatable           :: grid_xy(:,:,:)
    integer                     :: Npts, Nxpts, Nypts, idx, jdx, io
    logical                     :: exists
    
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


function solve_tridiag(a,b,c,f,J) result(u) !not necc most efficient b/c need to store abc, but for simplicity functionalized
!Solve trilinear eqn of the form Au = f, A is NxN, u and f are Nx1
!A is of form [ a_N    cN     0      0      0      ...    0  ]
!             [ b_Nm1  a_Nm1  c_Nm1  0      0      ...    0  ]
!             [ 0      b_Nm2  a_Nm2  c_Nm2  0      ...    0  ]
!               :      :      :      :      :             0  ]
!             [       0                     b_2    a_2    c_2]
!             [                             0      b_1    a_1]
!f is of form [f_N    f_Nm1   ... f_2    f_1]'            
!
!Slightly confusing formulation by mccormack. matrix index != index
!Performs inversion and back solve (Thomas algorithm?)
    implicit none
    integer, intent(in)                 :: J
    real(8), dimension(J),intent(in)    :: a,b,c,f
    real(8), dimension(J)               :: u
    
    integer                             :: jdx
    real(8)                             :: alj
    real(8),dimension(J+1)              :: gam,v!,u
    gam(J+1) = 0d0
    v(J+1) = 0d0
    
    !u(J+1) = 0d0
    do jdx = J,1,-1
        alj = a(jdx) - b(jdx)*gam(jdx+1)
        gam(jdx) = c(jdx)/alj
        v(jdx) = (f(jdx)-b(jdx)*v(jdx+1))/alj
    end do
    u(1) = v(1)
    do jdx = 2,J
        u(jdx) = v(jdx)-gam(jdx)*u(jdx-1)
    end do
end function solve_tridiag

function solve_mc(x,y,phi_IC,I,J,nstop,freestream) result(phi_n)
    implicit none
    integer, intent(in)                 :: I,J 
    real(8), dimension(I,J),intent(in)  :: x,y,phi_IC
    integer, intent(in)                 :: nstop
                                        
    real(8), dimension(I,J)             :: phi, A, mu
    real(8), dimension(J)               :: ad,bd,cd,fd
    real(8)                             :: dx1,dx2,dx3,dy1,dy2
    real(8)                             :: Rij, res
    integer                             :: idx,jdx,ndx
    
    !initalize phi, A,mu
    phi = phi_IC
    do idx = 2,I-1
        do jdx = 2,J-1
            A(idx,jdx) = 1d0-M_inf**2d0-(gam+1d0)*(M_inf**2)/V_inf * &
            (phi(idx+1,jdx)-phi(idx-1,jdx))/(x(idx+1,jdx)-x(idx-1,jdx))
            if (A(idx,jdx)>=0d0) then!faster way to do this than logicals? (needed?)
                mu(idx,jdx) = 0d0
            else
                mu(idx,jdx) = 1d0
            end if
        end do
    end do
    
    do ndx = 1,nstop
        do idx = 3,I-1 !TODO need special eqn at idx = 2, neglect supersonic term (idx-2)
            !assemble line
            do jdx = 2,J
                dx1 = x(idx+1,jdx) - x(idx,jdx)
                dx2 = x(idx,jdx)   - x(idx-1,jdx)
                dx3 = x(idx-1,jdx) - x(idx-2,jdx)
                dy1 = y(idx,jdx+1) - y(idx,jdx)
                dy2 = y(idx,jdx)   - y(idx,jdx-1)
                
                ad(jdx) = (2d0*(1d0-mu(idx,jdx))*A(idx,jdx))/(dx1+dx2) * (1d0/dx1+1d0/dx2)
                ad(jdx) = ad - (2d0*mu(idx-1,jdx)*A(idx-1,jdx)/(dx2+dx3) * 1d0/dx2
                ad(jdx) = ad + 2d0/(dy1+dy2)*(1d0/dy+1d0/dy2)
              
                bd(jdx) = -2d0/(dy1+dy2) * (1d0/dy1)
                cd(jdx) = -2d0/(dy1+dy2) * (1d0/dy2)
                
                fd(jdx) = (2d0*(1d0-mu(idx,jdx))*A(idx,jdx))/(dx1+dx2) * (phi(idx+1,jdx)/dx1+phi(idx-1,jdx)/dx2)
                fd(jdx) = fd(jdx) + (2d0*mu(idx-1,jdx)*A(idx-1,jdx)/(dx2+dx3) * &
                (-phi(idx-1,jdx)/dx2 - phi(idx-1,jdx)/dx3 +phi(idx-2,jdx)/dx3) 
            end do
            !boundary conditions
            ad(J) = 1d0
            bd(J) = 0d0 !DNE in real matrix
            cd(J) = 0d0
            fd(J) = phi_inf
            
            ad(1) = 1d0
            bd(1) = -1d0
            cd(1) = 0d0 !DNE in real mat
            if (x(idx,1)<0d0 .or. x(idx,1)>c) then  
                dphidy = 0d0
            else
                xb = r*cos(atan2((r-0.5e0*th),(x(idx,1)-c/2e0)))+0.5d0*c!body x (x is chord line)
                dphidy = M*sqrt(gama*P_inf/rho_inf) * (xb-c/2d0)/(r*sqrt(1d0-((xb**2d0-c/2d0)/(r))**2))
            end if
            f(1) = -dphidy * (y(idx,2)-y(idx,1))
            
            phi(idx,:) = solve_tridiag(ad,bd,cd,fd,J)
        end do
        !compute new A,mu (want consistent for resid, so done here)
        do idx = 2,I-1
            do jdx = 2,J-1
                A(idx,jdx) = 1d0-M_inf**2d0-(gam+1d0)*(M_inf**2)/V_inf * &
                (phi(idx+1,jdx)-phi(idx-1,jdx))/(x(idx+1,jdx)-x(idx-1,jdx))
                if (A(idx,jdx)>=0d0) then!faster way to do this than logicals? (needed?)
                    mu(idx,jdx) = 0d0
                else
                    mu(idx,jdx) = 1d0
                end if
            end do
        end do
        
        !compute and write residual
        do idx = 3,I-1 !TODO also need special eqn here for idx = 2
            do jdx = 2,J-1
                dx1 = x(idx+1,jdx) - x(idx,jdx)
                dx2 = x(idx,jdx)   - x(idx-1,jdx)
                dx3 = x(idx-1,jdx) - x(idx-2,jdx)
                dy1 = y(idx,jdx+1) - y(idx,jdx)
                dy2 = y(idx,jdx)   - y(idx,jdx-1)
                
                Rij =       (2d0*(1d0-mu(idx,jdx))*A(idx,jdx)) / (dx1+dx2) * &
                   ((phi(idx+1,jdx)-phi(idx,jdx))/dx1 - (phi(idx,jdx)-phi(idx-1,jdx))/dx2)
                Rij = Rij + (2d0*mu(idx,jdx)*A(idx,jdx)) / (dx2+dx3) * &
                   ((phi(idx,jdx+1)-phi(idx-1,jdx))/dx2 - (phi(idx-1,jdx)-phi(idx-2,jdx))/dx3)
                Rij = Rij + 2d0/(dy1+dy2) * &
                   ((phi(idx,jdx+1)-phi(idx,jdx))/dy1 - (phi(idx,jdx)-phi(idx,jdx-1))/dy2)
                if abs(Rij)> res then 
                    res = abs(Rij)
                end if
            end do
        end do
        print *,ndx,res
        
    end do

end function solve_mc 




program murman_cole
    implicit none
    
    integer :: io
    logical :: exists
    
    character(:), allocatable :: grid_path
    real(kdub) :: M
    real(kdub) :: curX
    integer :: nstop
    
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
    
    
end program murman_cole