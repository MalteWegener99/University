      
        subroutine hypot_normal(x, y, h)
        double precision :: x, y, h

        h = SQRT(x*x + y*y)
        end subroutine
      
        subroutine hypot_safer_fabs(x, y, h)
        double precision :: x, y, h

        if (x > y) then
            h = ABS(x) * SQRT(1.0 +(y/x)*(y/x))
        else if(x < y) then
            h = ABS(y) * SQRT(1.0 +(x/y)*(x/y))
        else
            h = SQRT(2.0) * ABS(x)
        end if
        end subroutine
