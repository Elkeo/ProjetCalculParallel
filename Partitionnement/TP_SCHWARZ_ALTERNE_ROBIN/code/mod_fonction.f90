module fonctions
use physical_parameter
implicit none

contains 
  function sol_exacte(a) result(b)
      implicit none
      double precision, intent(in) :: a
      double precision             :: b

      b = (a-x_min)*(a-x_max)*exp(-a)
      return
  end function


  function source_term(a) result(b)
      implicit none
      double precision, intent(in) :: a
      double precision             :: b

      b = -(a*a-5.d0*a+4.d0)*exp(-a)
      return
  end function

end module fonctions
