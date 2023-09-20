function kap_fun(kappa, D, JL, dy_min) result(f)
     implicit none
     integer, parameter :: kflt = SELECTED_REAL_KIND(9,10) !general float, 9 digits, max val 1e10

     real(kflt),intent(in) :: kappa
     real(kflt),intent(in) :: D
     real(kflt),intent(in) :: JL
     real(kflt),intent(in) :: dy_min

     real(kflt) :: f              

      f = dy_min - D*((exp(kappa*(1/(JL-1)))-1)/(exp(kappa)-1))
end function kap_fun

function dkap_fun(kappa, D, JL) result(fp)
     implicit none
     integer, parameter :: kflt = SELECTED_REAL_KIND(9,10) !general float, 9 digits, max val 1e10
     real(kflt),intent(in) :: kappa
     real(kflt),intent(in) :: D
     real(kflt),intent(in) :: JL
     real(kflt) :: fp
          
      fp = -D*(((-exp(kappa)*(exp(kappa*(1/(JL-1)))-1))/((exp(kappa)-1)**2)) &
              + ((exp(kappa*(1/(JL-1))))/((JL-1)*(exp(kappa)-1))))
end function dkap_fun


program murman_cole
      implicit none 
      
      integer, parameter :: kint = SELECTED_INT_KIND(16) !general int, 16 digits
      integer, parameter :: kflt = SELECTED_REAL_KIND(9,10) !general float, 9 digits, max val 1e10
      integer, parameter :: kdub = SELECTED_REAL_KIND(20,200) !general "double", 20 digits, max val 1e200
  
      !function return types
      real(kflt) :: kap_fun
      real(kflt) :: dkap_fun
      !control constants
      integer :: io
      
      !I/O files
      character(:), allocatable :: gridOutPath

      !parameters
      integer(kint), parameter :: grid_Npts = 51_kint !# grid points
      integer(kint), parameter :: kmax =10_kint; !num steps to conv kappa
      real(kflt), parameter :: c = 1.0_kflt !chord
      real(kflt), parameter :: th = 0.06_kflt !airfoil thickness (dfines dy_min)
      real(kflt), parameter :: ymin = 0.0_kflt !airfoil thickness (dfines dy_min)
      
      !real vars
      real(kflt) :: D !max y dimension
      real(kflt) :: dy_min !minimum y step
      real(kflt) :: kappa_k
      real(kflt) :: kappa_kp1
      real(kflt) :: kappa
      
      real(kflt), dimension(grid_Npts) :: grid_x
      real(kflt), dimension(grid_Npts) :: grid_y

      integer(kint) :: idx !general indexes. Should only use in loops
      integer(kint) :: jdx
      integer(kint) :: kdx
      
      !Assign vars explicitly
      gridOutPath = 'grid.txt'
      D = 50*c
      dy_min = th/10 
     
      print *, D, dy_min 
      !find kappa for grid spacing
      kappa_k = 1;
      do kdx = 1, kmax 
         kappa_kp1 = kappa_k -(kap_fun(kappa_k, D, real(grid_Npts,kflt), dy_min)/dkap_fun(kappa_k, D, real(grid_Npts,kflt)))
!         print *, "kappa=", kappa_k, ", f(kappa)=",kap_fun(kappa_k, D, real(grid_Npts,kflt), dy_min), &
!                 ", df(kappa)/dkappa=", dkap_fun(kappa_k, D, real(grid_Npts,kflt)) 
         kappa_k = kappa_kp1
      end do 
      kappa = kappa_kp1

      !Generate Grid with found kappa
      grid_y(1) = 0 
      do jdx=2,grid_Npts
          grid_y(jdx) = grid_y(jdx-1)+D*((exp(kappa*((real(jdx,kflt)-1)/(real(grid_Npts,kflt)-1)))-1)/(exp(kappa)-1))
          print *, grid_y(jdx) 
      end do 


!      open (newunit=io, file=gridOutPath, status="replace", action="write")
!      do idx=1, grid_Npts  
!          write(io, *) circ_xy(idx,1),  circ_xy(idx,2)
!      end do
      
!      close(io)
!      do idx=1, circ_Npts
          
!      end do



end program murman_cole
