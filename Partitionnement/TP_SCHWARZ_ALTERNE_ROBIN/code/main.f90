program fusion
   use physical_parameter
   use fonctions
   use BiCGStab_mod

   implicit none

   integer :: itschwz, maxschwarz, i
   integer :: alloc_stat
   real (kind=8), dimension(3) :: stencil !stencil du schema pour les conditions de bords Robin

   real (kind=8) :: error, errschwz, keep

   maxschwarz = 10000
   error      = 100.d0
   errschwz   = 1.e-8

   call read_input       ! Read input parameters
   allocate (X0(N0))
   allocate (X1(N1))
   allocate (U0(N0))
   allocate (U1(N1))
   allocate (RHS0(N0))
   allocate (RHS1(N1))
! Initialisation
   Do i =1, N0        !Premier sous-domaine Me == 0
      X0(i)   = i*h
      U0(i)   = 0.d0
   End DO
   Do i = N1, 1, -1    !Deuxième sous-domaine Me == 1
      X1(1+N1-i)   = x_max-i*h
      U1(1+N1-i)   = 0.d0
   End Do
   stencil=0.0d0

!Schwarz loop
   itschwz = 0
   Do while ((itschwz .le. maxschwarz) .and. (error .gt. errschwz))
      keep =U0(N0-1)
      RHS0 = second_membre(X0,stencil,0)
      U0   = BiCGStab(RHS0,0) ! On resoud le sous domaine 0
      !On met a jour le stencil avec la solution U0
      stencil(1) = U0(N0 - 2);stencil(2)=U0(N0 - 1);stencil(3)=U0(N0);
      !On resoud sous domaine 1
      RHS1 = second_membre(X1,stencil,1)
      U1   = BiCGStab(RHS1,1)
      !On met a jour le stencil avec la solution U1
      stencil(1) = U1(1) ;stencil(2) = U1(2);stencil(3)=U1(3);
      itschwz = itschwz + 1
      error = max(abs(U0(N0-1)-U1(1)), abs(U0(N0)-U1(2)))
   End Do

   write(*,*) "La convergence demande ", itschwz, " iterations de Schwarz"
   open(unit=10, file='Sol_000.dat', status='unknown')
   do i=1, N0
      write(10,*) X0(i), U0(i), sol_exacte(X0(i))
   End Do
   close(10)
   open(unit=11, file='Sol_001.dat', status='unknown')
   do i=1, N1
      write(11,*) X1(i), U1(i), sol_exacte(X1(i))
   End Do
   close(11)

   deallocate(X0,X1,U0,U1,RHS0,RHS1)

end program fusion
