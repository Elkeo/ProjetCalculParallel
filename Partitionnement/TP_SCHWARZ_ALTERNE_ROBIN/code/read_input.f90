subroutine read_input
use physical_parameter

implicit none

x_min = 0.d0
x_max = 1.d0

!Default parameters
! On veut 20 points sur le domaine 0, et 20 sur domaine 1, les extremites X_min, X_max ne sont pas inconnues Dirichlet classique
! Attention nous avons un recouvrement de 2 minimum
N = 38 
h = (x_max - x_min)/(N+1)

N0 = (N+2)/2
N1 = (N+2) - N0

open(unit=10, file='input_data', status='old')
read(10,*) alpha
read(10,*) beta
close(10)


write(*,*) 'N= ', N, ' h= ', h, ' N0= ', N0, ' N1= ', N1

end subroutine read_input
