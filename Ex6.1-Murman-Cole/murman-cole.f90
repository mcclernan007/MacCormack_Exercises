program murman_cole
      implicit none 
      
      integer, parameter :: kint = SELECTED_INT_KIND(16) !general int, 16 digits
      integer, parameter :: kflt = SELECTED_REAL_KIND(9,10) !general float, 9 digits, max val 1e10
      integer, parameter :: kdub = SELECTED_REAL_KIND(20,200) !general "double", 20 digits, max val 1e200

      real(kflt) :: pi
      
      real(kflt) :: airfoil_x
      real(kflt) :: airfoil_y
      real(kflt) :: airfoil_goalTh
      real(kflt) :: airfoil_th
      real(kflt) :: airfoil_chord

      integer(kint), parameter :: circ_Npts = 200_kint

      real(kflt), dimension(circ_Npts) :: circ_tht
      real(kflt) :: circ_dtht
      real(kflt), dimension(circ_Npts) :: circ_x
      real(kflt), dimension(circ_Npts) :: circ_y
      real(kflt), dimension(circ_Npts,2) :: circ_xy

      integer :: io

      integer(kint) :: idx !general indexes. Should only use in loops
      integer(kint) :: jdx
      integer(kint) :: kdx
    
      pi = 3.141592653589793238462643383279502_kflt


      airfoil_goalTh = 0.06
      airfoil_chord = 1.0

  
!generate unit quarter circle, shift by varied dy, scale st c=1, measure thickness. get close. Need array knowlege
      circ_dtht = (pi/2_kflt)/(real(circ_Npts)-1_kflt)
      circ_tht(1) = 0_kflt
      do idx=2_kint, circ_Npts
        circ_tht(idx) = circ_tht(idx-1)+circ_dtht
      end do 

      circ_x = 1*cos(circ_tht) 
      circ_y = 1*sin(circ_tht)
      circ_xy(:,1) = circ_x
      circ_xy(:,2) = circ_y
      
      open (newunit=io, file="out.txt", status="replace", action="write")
      do idx=1, circ_Npts  
          write(io, *) circ_xy(idx,1),  circ_xy(idx,2)
      end do
      
      close(io)
!      do idx=1, circ_Npts
          
!      end do
      


     !start by finding airfoil 


end program murman_cole
