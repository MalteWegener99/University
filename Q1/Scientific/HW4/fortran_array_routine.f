      
      subroutine multiply_array(a, b, c, n)
      
      integer a(10) , b(10) , c(10)
      integer :: n
      integer :: x
      
      do x = 1, n
        c(x) = a(x) * b(x)
      end do
      
      end subroutine
      
