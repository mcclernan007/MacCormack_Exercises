program murman_cole
      
      integer, parameter :: k_int = SELECTED_INT_KIND(16) !general int, 16 digits
      integer, parameter :: k_flt = SELECTED_REAL_KIND(9,10) !general float, 9 digits, max val 1e10
      integer, parameter :: k_dub = SELECTED_REAL_KIND(20,200) !general "double", 20 digits, max val 1e200

      real(k_flt) :: pi
      
      real(k_flt) :: airfoil_x
      real(k_flt) :: airfoil_y
      real(k_flt) :: airfoil_goalTh
      real(k_flt) :: airfoil_th
      real(k_flt) :: airfoil_chord

      real(k_flt) :: circ_x
      real(k_flt) :: circ_y
      integer(k_int) :: circ_Npts

      integer(k_int) :: idx !general indexes. Should only use in loops
      integer(k_int) :: jdx
      integer(k_int) :: kdx


      
      pi = 3.141592653589793238462643383279502_k_flt

      airfoil_goalTh = 0.06
      airfoil_chord = 1.0

      !generate unit circle, shift by varied dy, scale st c=1, measure thickness. get close. Need array knowlege
      circ_Npts = 200 
     
      circ_x = 
      do idx=1, circ_Npts
          
      end do
      


     !start by finding airfoil 


end program murman_cole
