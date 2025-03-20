module mod_Impl_PredVel
use mod_param
use mod_common
use mod_solver
implicit none
private
public Impl_PredVel
contains
subroutine Impl_PredVel(ustar,vstar,wstar)
real, intent(inout), dimension(0:,0:,0:) :: ustar,vstar,wstar
integer :: i,j,k,knum
real :: coef,alpha
real, dimension(1:kmax) :: A,B,C
real::nu_m

alpha = 0.5*dt
nu_m =0.5*(vis_1/rho_1+vis_2/rho_2)
coef = alpha*nu_m

do k=1,kmax
  A(k) = -coef*dzi*dzi
  C(k) = -coef*dzi*dzi
  B(k) = 1. + 2.*coef*dzi*dzi
enddo

B(1) = 1. + 3.*coef*dzi*dzi 
A(1) = 0.
B(kmax) = 1. + 3.*coef*dzi*dzi
C(kmax) = 0.
knum = kmax
call solver2d_Star(ustar,A,B,C,knum,coef)
call solver2d_Star(vstar,A,B,C,knum,coef)

B(1) = 1. + 2.*coef*dzi*dzi
A(1) = 0.
B(kmax) = 1. + 2.*coef*dzi*dzi
C(kmax) = 0.
C(kmax-1) = 0.
knum = kmax - 1
call solver2d_Star(wstar,A,B,C,knum,coef)

!
return
end subroutine Impl_PredVel
!
end module mod_Impl_PredVel

