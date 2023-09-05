program variables
      implicit none
        
      !declare vars with types
      integer :: amount
      real :: pi
      complex :: frequency
      character :: initial
      logical :: isOkay
        
      !assign variables
      amount = 10
      pi = 3.14
      frequency = (1,-0.5)
      initial = 'a'
      isOkay = .false.

      !writing to std out
      print *, 'the value of amount(integer) is: ', amount
      print *, 'thevalue of pi(real) is: ', pi
      print *, 'freq (complex) is :', frequency
      print *, 'initial(char) is:', initial
      print *, 'isOK(bool) is:', isOkay

      !stdin 


end program variables




