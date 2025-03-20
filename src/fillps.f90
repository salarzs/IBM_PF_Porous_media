module mod_fillps
use mod_param
use mod_common
use mod_common_IBM
implicit none
private
public fillps
contains
subroutine fillps(p)
real, intent(out), dimension(0:,0:,0:) :: p
real :: dti,dtidxi,dtidyi,dtidzi,rho0
integer :: i,j,k,im,jm,km
!
rho0=MIN(rho_1,rho_2)
dti = 1./dt
dtidxi = dti*dxi
dtidyi = dti*dyi
dtidzi = dti*dzi
do k=1,kmax
  km = k-1
  do j=1,jmax
    jm = j-1
    do i=1,imax
      p(i,j,k) = 0.
       im = i-1
       p(i,j,k) = ( &
                   ( dwdt(i,j,k)-dwdt(i,j,km))*dtidzi+ &
                   ( dvdt(i,j,k)-dvdt(i,jm,k))*dtidyi+ &
                   ( dudt(i,j,k)-dudt(im,j,k))*dtidxi )
 
       p(i,j,k) = p(i,j,k)*rho0 + (&
                ((1.0-rho0/(rhol(i,j,k)+rhol(i+1,j,k))*2.0)*(phat(i+1,j,k)-phat(i,j,k))*dxi- &
                 (1.0-rho0/(rhol(i,j,k)+rhol(i-1,j,k))*2.0)*(phat(i,j,k)-phat(i-1,j,k))*dxi)*dxi+ &
                ((1.0-rho0/(rhol(i,j,k)+rhol(i,j+1,k))*2.0)*(phat(i,j+1,k)-phat(i,j,k))*dyi- &
                 (1.0-rho0/(rhol(i,j,k)+rhol(i,j-1,k))*2.0)*(phat(i,j,k)-phat(i,j-1,k))*dyi)*dyi+ &
                ((1.0-rho0/(rhol(i,j,k)+rhol(i,j,k+1))*2.0)*(phat(i,j,k+1)-phat(i,j,k))*dzi- &
                 (1.0-rho0/(rhol(i,j,k)+rhol(i,j,k-1))*2.0)*(phat(i,j,k)-phat(i,j,k-1))*dzi)*dzi) 
   enddo
  enddo
enddo
!
return
end subroutine fillps
!
end module mod_fillps
