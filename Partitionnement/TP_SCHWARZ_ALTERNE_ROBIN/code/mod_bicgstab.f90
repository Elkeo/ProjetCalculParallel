!---------------------------------------------------------------------------------------!
!       REFERENCE :                                                                     !
!       The Improved BiCGStab Method for Large and Sparse Unsymmetric Linear            !
!       Systems on Parallel Distributed Memory Architectures; Yang et. al.              !
!                                                                                       !
!       NOTE:                                                                           !
!       This is the classic version of BICGStab. Stopping criteria of these code        !
!       needs additional treatment to avoid/detect failure condition!                   !
!       This code is tested with gFortran 5.01                                          !
!---------------------------------------------------------------------------------------!

module BiCGStab_mod
   implicit none

contains

   function BiCGStab(b, Me) result(x)
      implicit none

!--------------------------PARAMETER AND VARIABLE-------------------------------!
      real(kind=8), intent(in)                   :: b(:)
      real(kind=8), dimension(1:size(b, dim=1))   :: x
      integer, intent(in)                             :: Me

      real(kind=8), dimension(1:size(b, dim=1))   :: r, rs, v, p, s, t
      real(kind=8), parameter                     :: e = 1d-10
      real(kind=8)                                :: rho, rho_prev
      real(kind=8)                                :: alpha, omega, beta
      real(kind=8)                                :: norm_r, norm_b
      real(kind=8)                                :: summesion, temp

      integer                                         :: it, err, itmax = 100000
!------------------------END PARAMETER AND VARIABLE-----------------------------!

      !if(size(A, dim=1) /= size(A, dim=2)) stop &
      !"Error: Improper dimension of matrix A in BiCGStab."
      it = 0

!-------------------------------------------------------!
      x = 0.0d0                                      !-------> INITIAL GUESS
!-------------------------------------------------------!
      r = b - mon_matmul(x, Me)                            !-------> LINE 1
      rs = r                                          !
!-------------------------------------------------------!
      rho = 1.0d0; alpha = 1.0d0; omega = 1.0d0     !-------> LINE 2
!-------------------------------------------------------!
      v = 0.0d0; p = 0.0d0                          !-------> LINE 3
!                                                       !
      norm_r = dsqrt(dot_product(r, r))                !
      norm_b = dsqrt(dot_product(b, b))                !
!-------------------------------------------------------!

      do while ((norm_r .GT. e*norm_b) .AND. (it < itmax))        !-------> START OF LOOP

         !-------------------------------------------------------!
         rho_prev = rho                                      !-------> LINE 5
         rho = dot_product(rs, r)                        !
         !-------------------------------------------------------!
         beta = (rho/rho_prev)*(alpha/omega)           !-------> LINE 6
         !-------------------------------------------------------!
         p = r + beta*(p - omega*v)                 !-------> LINE 7
         !-------------------------------------------------------!
         v = mon_matmul(p, Me)                          !-------> LINE 8
         !-------------------------------------------------------!
         alpha = rho/dot_product(rs, v)                    !-------> LINE 9
         !-------------------------------------------------------!
         s = r - alpha*v                              !-------> LINE 10
         !-------------------------------------------------------!
         t = mon_matmul(s, Me)                          !-------> LINE 11
         !-------------------------------------------------------!
         omega = dot_product(t, s)/dot_product(t, t)        !-------> LINE 12
         !-------------------------------------------------------!
         x = x + alpha*p + omega*s                    !-------> LINE 13
         !-------------------------------------------------------!
         r = s - omega*t                              !-------> LINE 17
         !-------------------------------------------------------!
         norm_r = dsqrt(dot_product(r, r))
         norm_b = dsqrt(dot_product(b, b))

         it = it + 1

      end do                                                      !-------> END OF LOOP

      !print*, "Iteration Required :", it

      return
   end function BiCGStab

   function mon_matmul(b, Me) result(W)
      use physical_parameter

      implicit none
      real(kind=8), intent(in) :: b(:)
      integer, intent(in)       :: Me
      real(kind=8), dimension(1:size(b, dim=1)) :: W

      real(kind=8)             :: A11, A12, Aii, Aij
      integer                   :: i

      Aii = 2.d0/(h*h); Aij = -1.d0/(h*h); !Schema classique standard
      if (alpha .ne. 0) then
         A11 = Aii + 2*beta/(h*alpha); !Termes particuliers Dirichlet - pas necessairement sur la premiere ligne
         A12 = -Aii;
      else
         A12 = 0.0d0; !Termes particuliers Robin - pas necessairement sur la premiere ligne
         A11 = -Aij;
      end if

      if (Me .eq. 0) then
         W(1) = Aii*b(1) + Aij*b(2)                                        !Premiere ligne
         W(size(b, dim=1)) = A11*b(size(b, dim=1)) + A12*b(size(b, dim=1) - 1)  !Derniere ligne
      else
         W(1) = A11*b(1) + A12*b(2)                                        !Premiere ligne
         W(size(b, dim=1)) = Aii*b(size(b, dim=1)) + Aij*b(size(b, dim=1) - 1)  !Derniere ligne
      end if
      Do i = 2, size(b, dim=1) - 1
         W(i) = Aij*b(i - 1) + Aii*b(i) + Aij*b(i + 1)    !ligne standard
      End Do

      return
   end function mon_matmul

   function second_membre(xcoord, stencil, Me) result(RHS)
      use physical_parameter
      use fonctions
      implicit none
      integer, intent(in)                            :: Me
      real(kind=8), intent(in)                      :: xcoord(:)
      real(kind=8), dimension(3), intent(in)        :: stencil
      real(kind=8), dimension(1:size(xcoord, dim=1)) :: RHS

      integer :: i

      Do i = 1, size(xcoord, dim=1)
         RHS(i) = source_term(xcoord(i))
      End Do
      if (Me .eq. 0) then
         if (alpha .eq. 0.0d0) then
            RHS(size(xcoord, dim=1)) = stencil(2)/(h*h)
         else
            RHS(size(xcoord, dim=1)) = RHS(size(xcoord, dim=1)) + (stencil(3) - stencil(1))/(h*h) + 2.0d0*beta*stencil(2)/(h*alpha)
         end if
      else
         if (alpha .eq. 0.0d0) then
            RHS(1) = stencil(2)/(h*h)
         else
            RHS(1) = RHS(1) + (stencil(1) - stencil(3))/(h*h) + 2*beta*stencil(2)/(h*alpha)
         end if
      end if
      return
   end function second_membre

end module BiCGStab_mod
