function kap_fun(kappa, D, JL, dy_min) result(f)
     implicit none
     integer, parameter :: kflt = SELECTED_REAL_KIND(9,10) !general float, 9 digits, max val 1e10
     integer, parameter :: kint = SELECTED_INT_KIND(16) !general int, 16 digits

     real(kflt),intent(in) :: kappa
     real(kflt),intent(in) :: D
     integer(kint),intent(in) :: JL
     real(kflt),intent(in) :: dy_min

     real(kflt) :: f              

      f = dy_min - D*(exp(kappa*(1/(JL-1))-1)/(exp(kappa)-1))
end function kap_fun

function dkap_fun(kappa, D, JL) result(fp)
     implicit none
     integer, parameter :: kflt = SELECTED_REAL_KIND(9,10) !general float, 9 digits, max val 1e10
     integer, parameter :: kint = SELECTED_INT_KIND(16) !general int, 16 digits

     real(kflt),intent(in) :: kappa
     real(kflt),intent(in) :: D
     integer(kint),intent(in) :: JL
     real(kflt) :: fp
          
      fp = -D*((-exp(kappa)*(exp(kappa*(1/(JL-1)))-1))/((exp(kappa)-1)**2)&
              + (exp(kappa*(1/(JL-1))))/((JL-1)*(exp(kappa)-1)))
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
      integer(kint), parameter :: kmax = 50; !num steps to conv kappa
      real(kflt), parameter :: c = 1.0_kflt !chord
      real(kflt), parameter :: th = 0.06 !airfoil thickness (dfines dy_min)
      
      !real vars
      real(kflt) :: D !max y dimension
      real(kflt) :: dy_min !minimum y step
      real(kflt) :: kappa_prv
      real(kflt) :: kappa_cur
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
      
      !find kappa for grid spacing
      kappa_prv = 1;
      do idx = 1, kmax 
         kappa_cur = kappa_prv-kap_fun(kappa_prv, D, grid_Npts, dy_min)/dkap_fun(kappa_prv, D, grid_Npts)
         print *, kappa_prv
         kappa_prv = kappa_cur
      end do 
      kappa = kappa_cur

      !Generate Grid with spacing 
      
     
!      open (newunit=io, file=gridOutPath, status="replace", action="write")
!      do idx=1, grid_Npts  
!          write(io, *) circ_xy(idx,1),  circ_xy(idx,2)
!      end do
      
!      close(io)
!      do idx=1, circ_Npts
          
!      end do



end program murman_cole
