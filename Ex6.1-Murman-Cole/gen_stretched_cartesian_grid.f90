function kap_fun(kappa, D, JL, dy_min) result(f)
     implicit none
	 integer, parameter :: kint = SELECTED_INT_KIND(16) !general int, 16 digits
     integer, parameter :: kflt = SELECTED_REAL_KIND(9,10) !general float, 9 digits, max val 1e10

     real(kflt),intent(in) :: kappa
     real(kflt),intent(in) :: D
     integer(kint),intent(in) :: JL
     real(kflt),intent(in) :: dy_min

     real(kflt) :: f              

      f = dy_min - D*((exp(kappa*(1/(real(JL,kflt)-1)))-1)/(exp(kappa)-1))
end function kap_fun

function dkap_fun(kappa, D, JL) result(fp)
     implicit none
	 integer, parameter :: kint = SELECTED_INT_KIND(16) !general int, 16 digits
     integer, parameter :: kflt = SELECTED_REAL_KIND(9,10) !general float, 9 digits, max val 1e10
     real(kflt),intent(in) :: kappa
     real(kflt),intent(in) :: D
     integer(kint),intent(in) :: JL
     real(kflt) :: fp
          
      fp = -D*(((-exp(kappa)*(exp(kappa*(1/(real(JL,kflt)-1)))-1))/((exp(kappa)-1)**2)) &
              + ((exp(kappa*(1/(real(JL,kflt)-1))))/((JL-1)*(exp(kappa)-1))))
end function dkap_fun

function find_kappa(D, JL, d_min) result (kappa)
      implicit none
      integer, parameter :: kint = SELECTED_INT_KIND(16) !general int, 16 digits
	  integer, parameter :: kflt = SELECTED_REAL_KIND(9,10) !general float, 9 digits, max val 1e10
	  
      integer(kint), parameter :: kmax =10_kint; !num steps to conv kappa
	  
      real(kflt),intent(in) :: D
	  integer(kint),intent(in) :: JL
	  real(kflt),intent(in) :: d_min
	  real(kflt) :: kappa
	  
	  real(kflt) :: kappa_k
      real(kflt) :: kappa_kp1
	  integer(kint) :: kdx
	  
	  !function return types
	  real(kflt) :: kap_fun
      real(kflt) :: dkap_fun
	  
	  
      kappa_k = 1_kflt;
      do kdx = 1, kmax 
         kappa_kp1 = kappa_k - (kap_fun(kappa_k, D, JL, d_min)/dkap_fun(kappa_k, D,JL))
!         print *, "kappa=", kappa_k, ", f(kappa)=",kap_fun(kappa_k, D, real(grid_Npts,kflt), dy_min), &
!                 ", df(kappa)/dkappa=", dkap_fun(kappa_k, D, real(grid_Npts,kflt)) 
         kappa_k = kappa_kp1
      end do 
      kappa= kappa_kp1
end function find_kappa

program gen_grid
      implicit none 
      
      integer, parameter :: kint = SELECTED_INT_KIND(16) !general int, 16 digits
      integer, parameter :: kflt = SELECTED_REAL_KIND(9,10) !general float, 9 digits, max val 1e10
      integer, parameter :: kdub = SELECTED_REAL_KIND(20,200) !general "double", 20 digits, max val 1e200
  
      !function return types
      
	  real(kflt) :: find_kappa
      
	  !control constants
      integer :: io
      
      !I/O files
      character(:), allocatable :: gridOutPath

      !parameters
      integer(kint), parameter :: grid_Npts = 51_kint !# grid points
	  integer(kint), parameter :: grid_airfoilNpts = 21_kint !# grid const dx points on af
	  integer(kint), parameter :: grid_xStrechNpts = (grid_Npts-grid_airfoilNpts)/2 !#grid pts to stretch 
	  
      real(kflt), parameter :: c = 1.0_kflt !chord
      real(kflt), parameter :: th = 0.06_kflt !airfoil thickness (dfines dy_min)
      real(kflt), parameter :: ymin = 0.0_kflt !airfoil thickness (dfines dy_min)
      
      !real vars
      real(kflt) :: D !max y dimension
	  real(kflt) :: dx_min !min x step 
      real(kflt) :: dy_min !min y step

      real(kflt) :: kappa_x !final val kappa for x dir
	  real(kflt) :: kappa_y !final val kappa for y dir
	  
      real(kflt), dimension(grid_xStrechNpts) :: grid_x_str
      real(kflt), dimension(grid_Npts) :: grid_x
      real(kflt), dimension(grid_Npts) :: grid_y

      integer(kint) :: idx !general indexes. Should only use in loops
      integer(kint) :: jdx
      integer(kint) :: kdx
      
      !Assign vars explicitly
      gridOutPath = 'grid.txt'
      D = 50*c
	  dx_min = c/real(grid_airfoilNpts,kflt)
      dy_min = th/10 
     
      print *, D, dx_min, dy_min
      !find kappa for grid spacing
      kappa_x = find_kappa(D, grid_xStrechNpts, dx_min)
	  kappa_y = find_kappa(D, grid_Npts,dy_min)

      !Generate Grid with found kappa
	  !y
      grid_y(1) = 0_kflt
      do jdx=2,grid_Npts
          grid_y(jdx) = grid_y(1)+D*((exp(kappa_y*((real(jdx,kflt)-1)/(real(grid_Npts,kflt)-1)))-1)/(exp(kappa_y)-1))
!          print *, grid_y(jdx) 
      end do 
	  
	  !x
	  !constant spacing portion
	  grid_x(grid_xStrechNpts+1) = 0_kflt
	  do idx= grid_xStrechNpts+2, grid_xStrechNpts+grid_airfoilNpts, 1 
			grid_x(idx) = grid_x(idx-1) + c/real(grid_airfoilNpts-1,kflt)
	  end do
	  !stretched portion right
	  grid_x_str(1) = 0_kflt+dx_min;
	  
	  do idx=2, grid_xStrechNpts
         grid_x_str(idx) = grid_x_str(1)+D*((exp(kappa_x*((real(idx,kflt)-1)/ &
		 (real(grid_xStrechNpts,kflt)-1)))-1)/(exp(kappa_x)-1))
        ! print *, grid_x_str(idx) 
      end do
	  
	  grid_x(1:grid_xStrechNpts) = -1*grid_x_str(grid_xStrechNpts:1:-1)
	  grid_x(grid_xStrechNpts+grid_airfoilNpts+1:grid_Npts) = c + grid_x_str
	  
	  open (newunit=io, file=gridOutPath, status="replace", action="write")
	  write(io,*) grid_Npts*grid_Npts
	  do idx=1, grid_Npts
		 do jdx=1, grid_Npts
		 	write(io, *) grid_x(idx), " ", grid_y(jdx)
		 end do
	   end do
      
!      do idx=1, circ_Npts
          
!      end do



end program gen_grid