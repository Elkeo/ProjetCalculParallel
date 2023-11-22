module physical_parameter

!Numerical domain
double precision :: x_min, x_max

!Numerical parameters
integer          :: N !Total number of unknowns
integer          :: N0, N1 ! Number of unknowns for me = 0, me =1
double precision :: h
double precision :: alpha, beta

double precision, dimension(:), allocatable :: X0,X1,U0,U1,RHS0,RHS1


end module physical_parameter
