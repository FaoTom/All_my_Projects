!Week 1 - Exercise 2, Faorlin Tommaso 

program number_precision

    implicit none

    !variables for (a)
    integer*2 num_1_short
    integer*2 num_2_short
    integer*4 num_1_int
    integer*4 num_2_int

    !variables for (b)
    real*4 pi32_singlep
    real*4 var2_singlep
    real*8 pi32_doublep
    real*8 var2_doublep

    !definitions for (a)
    num_1_short = 2000000 !-31616 perchè prende i primi 16 bit dall'LSB
    num_2_short = 1
    num_1_int  = 2000000
    num_2_int  = 1

    !definitions for (b)
    pi32_singlep = 2.0*ASIN(1.0)*1e32 !or equivalently 4.0*ATAN(1.0)
    var2_singlep = SQRT(2.0)*1e21
    pi32_doublep = 2.D0*DASIN(1.D0)*1e32 !or equivalently 4.D0*DATAN(1.D0)
    var2_doublep = DSQRT(2.D0)*1e21
    
    print *, "(a) 2e6+1 with integer*2: " , num_1_short + num_2_short
    print *, "(a) 2e6+1 with integer*4: " , num_1_int + num_2_int
    print *, "(b) PIe32 + sqrt(2)e21 in single precision: " , pi32_singlep + var2_singlep
    print *, "(b) PIe32 + sqrt(2)e21 in double precision: " , pi32_doublep + var2_doublep

end program number_precision