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

subroutine get_freestream(M_inf,V_inf,rho_inf,gam,phi_inf,a_inf,P_inf)
    !nongeneral compact way to do this. Would want real input method 
    implicit none
    real(8), intent(in)  :: M_inf
    real(8), intent(out) :: V_inf, rho_inf, gam, phi_inf
    real(8), intent(out) :: a_inf, P_inf
    
    !assumed nondim
    V_inf = 1d0
    rho_inf = 1d0
    gam = 1.4d0 !assumed air, isent
    phi_inf = 1d0 !phihat_inf
    !calculated
    a_inf =V_inf/M_inf
    P_inf = 1d0/(gam*M_inf**2d0)
end subroutine get_freestream

subroutine get_airfoil(c,th,r)
    !nongeneral compact way to do this. Would want real input method 
    real(8),intent(out) :: c,th,r
    c = 1d0
    th = 0.06d0
    r = (c**2d0+th**2d0)/(4d0*th)
end subroutine get_airfoil

function get_IC(x,y,M_inf,I,J) result(phi_IC)
    implicit none
    integer, intent(in)                 :: I,J 
    real(8), dimension(I,J),intent(in)  :: x,y
    real(8), intent(in)                 :: M_inf
    real(8), dimension(I,J)             :: phi_IC
    
    real(8)                             :: V_inf,rho_inf,gam,phi_inf,a_inf,P_inf
    real(8)                             :: c,th,r, xb
    real(8)                             :: dphidy
    integer                             :: idx
    
    call get_freestream(M_inf,V_inf,rho_inf,gam,phi_inf,a_inf,P_inf)
    call get_airfoil(c,th,r)
    
    phi_IC(:,:) = phi_inf    
    do idx = 1,I
        if (x(idx,1)<0d0 .or. x(idx,1)>c) then  
            dphidy = 0d0 !always this, but could be given invalid grid
        else 
            !xb = r*cos(atan2((r-0.5d0*th),(x(idx,1)-c/2d0)))+0.5d0*c!body x (x is chord line)
            xb = x(idx,1)
            !dphidy = V_inf * (xb-c/2d0)/(r*sqrt(1d0-(((xb-c/2d0)**2d0)/(r**2d0))))
            dphidy = V_inf * -1d0*((xb-0.5d0*c)/(sqrt(r**2d0 - 0.25d0*(c-2d0*xb)**2)))
        end if
        phi_IC(idx,1) = phi_IC(idx,2)-dphidy * (y(idx,2)-y(idx,1))
    end do
end function get_IC

function solve_mc(x,y,phi_IC,I,J,nstop,M_inf,resid_path) result(phi)
    implicit none
    integer, intent(in)                 :: I,J 
    real(8), dimension(I,J),intent(in)  :: x,y,phi_IC
    real(8), intent(in)                 :: M_inf
    integer, intent(in)                 :: nstop
    character(*), intent(in)            :: resid_path
    real(8), dimension(I,J)             :: phi

    real(8), dimension(J)               :: ad,bd,cd,fd
    real(8)                             :: A_ij, A_im1j, mu_ij, mu_im1j
    real(8)                             :: V_inf,rho_inf,gam,phi_inf,a_inf,P_inf
    real(8)                             :: c,th,r, xb
    real(8)                             :: dx1,dx2,dx3,dy1,dy2
    real(8)                             :: dphidy
    real(8)                             :: Rij, res
    integer                             :: idx,jdx,ndx
    integer                             :: io
    
    interface
        function solve_tridiag(a,b,c,f,J) result(u) 
            integer, intent(in)             :: J
            real(8), dimension(J),intent(in)   :: a,b,c,f
            real(8), dimension(J)              :: u
        end function solve_tridiag
    end interface
    
    resmap(:,:) = 0d0
    open (newunit=io, file=resid_path, status="replace", action="write")
    
    
    call get_freestream(M_inf,V_inf,rho_inf,gam,phi_inf,a_inf,P_inf)
    call get_airfoil(c,th,r)
    
    
    phi = phi_IC
    
    do ndx = 1,nstop
        phi_prev = phi
        idx = 2 !special eqn at idx=2, neglect supersonic term (any w/mu(idx-1))
        do jdx = 2,J-1 !assemble line
            dx1 = x(idx+1,jdx) - x(idx,jdx)
            dx2 = x(idx,jdx)   - x(idx-1,jdx)
            dy1 = y(idx,jdx+1) - y(idx,jdx)
            dy2 = y(idx,jdx)   - y(idx,jdx-1)
            
            A_ij = 1d0-M_inf**2d0-(gam+1d0)*(M_inf**2)/V_inf * &
            (phi(idx+1,jdx)-phi(idx-1,jdx))/(x(idx+1,jdx)-x(idx-1,jdx))
            if (A_ij >= 0d0) then
                mu_ij = 0d0
            else
                mu_ij = 1d0
            end if
            
            ad(jdx) = (2d0*(1d0-mu_ij)*A_ij)/(dx1+dx2) * (1d0/dx1+1d0/dx2)
            ad(jdx) = ad(jdx) + 2d0/(dy1+dy2)*(1d0/dy1+1d0/dy2)
          
            bd(jdx) = -2d0/(dy1+dy2) * (1d0/dy1)
            cd(jdx) = -2d0/(dy1+dy2) * (1d0/dy2)
            
            fd(jdx) = (2d0*(1d0-mu_ij)*A_ij)/(dx1+dx2) * (phi(idx+1,jdx)/dx1+phi(idx-1,jdx)/dx2)
        end do
        !boundary conditions
        ad(J) = 1d0
        bd(J) = 0d0 !DNE in real matrix
        cd(J) = 0d0
        fd(J) = phi_inf
        
        ad(1) = -1d0
        bd(1) = 1d0
        cd(1) = 0d0 !DNE in real mat
        if (x(idx,1)<0d0 .or. x(idx,1)>c) then  
            dphidy = 0d0 !always this, but could be given invalid grid
        else 
            !xb = r*cos(atan2((r-0.5e0*th),(x(idx,1)-c/2e0)))+0.5d0*c!body x (x is chord line)
            xb = x(idx,1)
            dphidy = V_inf * -1d0*((xb-0.5d0*c)/(sqrt(r**2d0 - 0.25d0*(c-2d0*xb)**2)))
        end if
        fd(1) = dphidy * (y(idx,2)-y(idx,1))
        
        phi(idx,:) = solve_tridiag(ad,bd,cd,fd,J)
        do idx = 3,I-1 !rest of idx, w/full formulation
            !assemble line
            do jdx = 2,J-1
                dx1 = x(idx+1,jdx) - x(idx,jdx)
                dx2 = x(idx,jdx)   - x(idx-1,jdx)
                dx3 = x(idx-1,jdx) - x(idx-2,jdx)
                dy1 = y(idx,jdx+1) - y(idx,jdx)
                dy2 = y(idx,jdx)   - y(idx,jdx-1)
                
                
                A_ij = 1d0-M_inf**2d0 - (gam+1d0)*(M_inf**2)/V_inf * &
                (phi(idx+1,jdx)-phi(idx-1,jdx))/(x(idx+1,jdx)-x(idx-1,jdx))
                if (A_ij >= 0d0) then
                    mu_ij = 0d0
                else
                    mu_ij = 1d0
                end if
                A_im1j = 1d0-M_inf**2d0-(gam+1d0)*(M_inf**2)/V_inf * &
                (phi(idx,jdx)-phi(idx-2,jdx))/(x(idx,jdx)-x(idx-2,jdx))
                if (A_im1j >= 0d0) then
                    mu_im1j = 0d0
                else
                    mu_im1j = 1d0
                end if
                
                ad(jdx) = (2d0*(1d0-mu_ij)*A_ij)/(dx1+dx2) * (1d0/dx1+1d0/dx2)
                ad(jdx) = ad(jdx) - ((2d0*mu_im1j*A_im1j)/(dx2+dx3)) * (1d0/dx2)
                ad(jdx) = ad(jdx) + (2d0/(dy1+dy2)) * (1d0/dy1+1d0/dy2)
              
                bd(jdx) = -2d0/(dy1+dy2) * (1d0/dy1)
                cd(jdx) = -2d0/(dy1+dy2) * (1d0/dy2)

                fd(jdx) = (2d0*(1d0-mu_ij)*A_ij)/(dx1+dx2) * (phi(idx+1,jdx)/dx1+phi(idx-1,jdx)/dx2)
                fd(jdx) = fd(jdx) + (2d0*mu_im1j*A_im1j)/(dx2+dx3) * &
                (-phi(idx-1,jdx)/dx2 - phi(idx-1,jdx)/dx3 +phi(idx-2,jdx)/dx3) 

            end do
            !boundary conditions
            ad(J) = 1d0
            bd(J) = 0d0 !DNE in real matrix
            cd(J) = 0d0
            fd(J) = phi_inf
            
            ad(1) = -1d0
            bd(1) = 1d0
            cd(1) = 0d0 !DNE in real mat
            if (x(idx,1)<0d0 .or. x(idx,1)>c) then  
                dphidy = 0d0
            else
                !xb = r*cos(atan2((r-0.5d0*th),(x(idx,1)-c/2d0)))+0.5d0*c!body x (x is chord line)
                xb = x(idx,1)
                dphidy = V_inf * -1d0*((xb-0.5d0*c)/(sqrt(r**2d0 - 0.25d0*(c-2d0*xb)**2)))
            end if
            fd(1) = dphidy * (y(idx,2)-y(idx,1))
            
            phi(idx,:) = solve_tridiag(ad,bd,cd,fd,J)
        end do
        
        !compute and write residual
        res = -1d0
        idx = 2
        do jdx = 2,J-1 !special eqn at idx=2, neglect supersonic term (any w/mu(idx-1))
            A_ij = 1d0-M_inf**2d0-(gam+1d0)*(M_inf**2)/V_inf * &
            (phi(idx+1,jdx)-phi(idx-1,jdx))/(x(idx+1,jdx)-x(idx-1,jdx))
            if (A_ij >= 0d0) then
                mu_ij = 0d0
            else
                mu_ij = 1d0
            end if
                       
            Rij =       ((2d0*(1d0-mu_ij)*A_ij) / (dx1+dx2)) * &
                       ((phi(idx+1,jdx)-phi(idx,jdx))/dx1 - (phi(idx,jdx)-phi(idx-1,jdx))/dx2)
            Rij = Rij + (2d0/(dy1+dy2)) * &
                       ((phi(idx,jdx+1)-phi(idx,jdx))/dy1 - (phi(idx,jdx)-phi(idx,jdx-1))/dy2)
            resmap(idx,jdx) = Rij
            if (abs(Rij)> res) then 
                res = abs(Rij)
                idx_resmax = idx
                jdx_resmax = jdx
            end if
        end do
        do idx = 3,I-1 !rest of idx, w/full formulation
            do jdx = 2,J-1
                dx1 = x(idx+1,jdx) - x(idx,jdx)
                dx2 = x(idx,jdx)   - x(idx-1,jdx)
                dx3 = x(idx-1,jdx) - x(idx-2,jdx)
                dy1 = y(idx,jdx+1) - y(idx,jdx)
                dy2 = y(idx,jdx)   - y(idx,jdx-1)
                
                A_ij = 1d0-M_inf**2d0-(gam+1d0)*(M_inf**2)/V_inf * &
                (phi(idx+1,jdx)-phi(idx-1,jdx))/(x(idx+1,jdx)-x(idx-1,jdx))
                if (A_ij >= 0d0) then
                    mu_ij = 0d0
                else
                    mu_ij = 1d0
                end if
                A_im1j = 1d0-M_inf**2d0-(gam+1d0)*(M_inf**2)/V_inf * &
                (phi(idx,jdx)-phi(idx-2,jdx))/(x(idx,jdx)-x(idx-2,jdx))
                if (A_im1j >= 0d0) then
                    mu_im1j = 0d0
                else
                    mu_im1j = 1d0
                end if
                
                Rij =       ((2d0*(1d0-mu_ij)*A_ij) / (dx1+dx2)) * &
                   ((phi(idx+1,jdx)-phi(idx,jdx))/dx1 - (phi(idx,jdx)-phi(idx-1,jdx))/dx2)
                Rij = Rij + ((2d0*mu_im1j*A_im1j) / (dx2+dx3)) * &
                   ((phi(idx,jdx)-phi(idx-1,jdx))/dx2 - (phi(idx-1,jdx)-phi(idx-2,jdx))/dx3)
                Rij = Rij + (2d0/(dy1+dy2)) * &
                   ((phi(idx,jdx+1)-phi(idx,jdx))/dy1 - (phi(idx,jdx)-phi(idx,jdx-1))/dy2)
                resmap(idx,jdx) = Rij
                if (abs(Rij)> res) then 
                    res = abs(Rij)
                end if
            end do
        end do
        print *,ndx,res
        write(io,*) ndx, res
    end do
    close(io)
    
    open (newunit=io_deb3, file="res_map.dat", status="replace", action="write")
    write(io_deb3,*) I*J,I,J
    do idx=1, I
        do jdx=1,J
            write(io_deb3, *) x(idx,jdx),"  ", y(idx,jdx),"  ", resmap(idx,jdx)
        end do
    end do
    close(io_deb3)
    
    
end function solve_mc 

!subroutine get_uv(x,y, I,J, phi, u,v)   
!    integer, intent(in)                     :: I,J
!    real(8), dimension(I,J), intent(in)     :: x,y,phi
!    real(8), dimension(I,J), intent(out)    :: u,v
!    
!    integer :: idx,jdx
!    !TODO probably cant do centered differencing everywhere. find A, mu again or something
!    do idx = 2,I-1
!        do jdx = 1,J
!            u(idx,jdx) = phi(idx+1,idx-1)
!        end do
!    end do
!end subroutine get_uv

subroutine output_result(path, x,y, phi,I,J)
    implicit none
    character(*), intent(in)            :: path
    integer, intent(in)                 :: I,J
    real(8), dimension(I,J), intent(in)    :: phi,x,y
    
    integer :: io
    integer :: idx,jdx
    
    open (newunit=io, file=path, status="replace", action="write")
    write(io,*) I*J,I,J
    do idx=1, I
        do jdx=1,J
            write(io, *) x(idx,jdx),"  ", y(idx,jdx),"  ", phi(idx,jdx)
        end do
    end do
    close(io)
end subroutine output_result


program murman_cole
    implicit none
    
    real(8), allocatable   :: xy(:,:,:)
    real(8), allocatable   :: x(:,:),y(:,:), phi(:,:), phi_IC(:,:),u(:,:),v(:,:)
    real(8)                :: V_inf,rho_inf,gam,phi_inf,a_inf,P_inf
    character(:), allocatable :: grid_path
    real(8)                   :: M_inf
    integer                   :: nstop,I,J
    
    interface 
        function read_grid_file(grid_path) result(grid_xy)
            character(*), intent(in)    :: grid_path
            real(8), allocatable           :: grid_xy(:,:,:)
        end function read_grid_file
        function get_IC(x,y,M_inf,I,J) result(phi_IC)
            integer, intent(in)                 :: I,J 
            real(8), dimension(I,J),intent(in)  :: x,y
            real(8), intent(in)                 :: M_inf
            real(8), dimension(I,J)             :: phi_IC
        end function get_IC
        function solve_mc(x,y,phi_IC,I,J,nstop,M_inf,resid_path) result(phi)
            integer, intent(in)                 :: I,J 
            real(8), dimension(I,J),intent(in)  :: x,y,phi_IC
            real(8), intent(in)                 :: M_inf
            integer, intent(in)                 :: nstop
            character(*), intent(in)            :: resid_path
            real(8), dimension(I,J)             :: phi
        end function solve_mc
        
    end interface
    
    !M = 0.908_kdub
    
    grid_path = "./grid.txt"
    
    xy = read_grid_file("grid.txt") 
    I = size(xy,2)
    J = size(xy,3)
    
    allocate(x(I,J))
    allocate(y(I,J))
    allocate(phi(I,J))
    allocate(phi_IC(I,J))
    
    x = xy(1,:,:)
    y = xy(2,:,:)
    
    !!Case 1 #fully subsonic
    M_inf = 0.735d0
    nstop = 400 
    !nstop = 32000 !for convergence
    phi_IC = get_IC(x,y,M_inf,I,J)
    call output_result("phiIC-C1-M0p735.dat", x,y, phi_IC,I,J)
    
    phi = solve_mc(x,y,phi_IC,I,J,nstop,M_inf,"resid-C1-M0p735.dat")
    call output_result("phi-C1-M0p735.dat", x,y, phi,I,J)
    
    !Case 2 #actually transonic
    M_inf = 0.908d0
    nstop = 400
    !nstop = 27000 !for convergence
    phi_IC = get_IC(x,y,M_inf,I,J)
    call output_result("phiIC-C2-M0p908.dat", x,y, phi_IC,I,J)
    
    phi = solve_mc(x,y,phi_IC,I,J,nstop,M_inf,"resid-C2-M0p908.dat")
    call output_result("phi-C2-M0p908.dat", x,y, phi,I,J)
    !call get_uv(x,y, I,J, phi, u,v)   
    !call output_cp("cp-0-initial_conditions.dat",xy,u,v,I,J,M,gama,P_inf,rho_inf)
    
    
    
    
end program murman_cole