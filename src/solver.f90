module mod_solver
use mod_param
use mod_common
implicit none
private
public solver2d,solver2d_Star
contains
subroutine solver2d(pz,helm)
use decomp_2d
implicit none
real, intent(inout), dimension(1:,1:,1:) :: pz
real,intent(in) :: helm
real, dimension(itot/dims(1),jtot,ktot/dims(2)) :: py
real, dimension(itot,jtot/dims(1),ktot/dims(2)) :: px
real :: bb
real :: z,d(imax,jmax,kmax)
real :: di(itot),dj(jtot)
integer :: i,j,k
do k=1,kmax
b(k) = b(k) + helm
enddo
call transpose_z_to_y(pz,py)
call transpose_y_to_x(py,px)
!
do k=1,xsize(3)
  do j=1,xsize(2)
    call vrfftf(1,itot,px(1:itot,j,k),di,1,wi)
  enddo
enddo
!
call transpose_x_to_y(px,py)
!
do k=1,ysize(3)
  do i=1,ysize(1)
    call vrfftf(1,jtot,py(i,1:jtot,k),dj,1,wj)
  enddo
enddo
!
call transpose_y_to_z(py,pz)
!
do j=1,jmax
  do i=1,imax
    z        = 1./(b(1)+xyrt(i,j))
    d(i,j,1) = c(1)*z
    pz(i,j,1) = pz(i,j,1)*z
  enddo
enddo
do k=2,kmax-1
   do j=1,jmax
     do i=1,imax
       bb       = b(k)+xyrt(i,j)
       z        = 1./(bb-a(k)*d(i,j,k-1))
       d(i,j,k) = c(k)*z
       pz(i,j,k) = (pz(i,j,k)-a(k)*pz(i,j,k-1))*z
     enddo
  enddo
enddo
do j=1,jmax
  do i=1,imax
    bb       = b(kmax)+xyrt(i,j)
    z        = bb-a(kmax)*d(i,j,kmax-1)
    if(z.ne.0.) then
      pz(i,j,kmax) = (pz(i,j,kmax)-a(kmax)*pz(i,j,kmax-1))/z
    else
      pz(i,j,kmax) =0.
    endif
  enddo
enddo

do k=kmax-1,1,-1
  do j=1,jmax
    do i=1,imax
      pz(i,j,k) = pz(i,j,k)-d(i,j,k)*pz(i,j,k+1)
    enddo
  enddo
enddo

call transpose_z_to_y(pz,py)
do k=1,ysize(3)
  do i=1,ysize(1)
    call vrfftb(1,jtot,py(i,1:jtot,k),dj,1,wj)
  enddo
enddo
!
call transpose_y_to_x(py,px)
!
do k=1,xsize(3)
  do j=1,xsize(2)
    call vrfftb(1,itot,px(1:itot,j,k),di,1,wi)
  enddo
enddo
!
call transpose_x_to_y(px,py)
call transpose_y_to_z(py,pz)
!
return
end subroutine solver2d

subroutine solver2d_Star(pz,aaa,bbb,ccc,kkmax,coef)
use decomp_2d
implicit none
real, intent(inout), dimension(1:,1:,1:) :: pz
real, intent(in), dimension(1:kmax) :: aaa,bbb,ccc
real, intent(in) :: coef
integer, intent(in) :: kkmax
real, dimension(itot/dims(1),jtot,ktot/dims(2)) :: py
real, dimension(itot,jtot/dims(1),ktot/dims(2)) :: px
real :: bb
real :: z,d(imax,jmax,kmax)
real :: di(itot),dj(jtot)
integer :: i,j,k
!
call transpose_z_to_y(pz,py)
call transpose_y_to_x(py,px)
!
!$omp parallel default(none) &
!$omp& shared(px,xsize) private(i,j,k,di) &
!$omp&firstprivate(wi)
!$omp do
do k=1,xsize(3)
  do j=1,xsize(2)
    call vrfftf(1,itot,px(1:itot,j,k),di,1,wi)
  enddo
enddo
!$omp end parallel
!
call transpose_x_to_y(px,py)
!
!$omp parallel default(none) &
!$omp& shared(py,ysize) private(i,j,k,dj) &
!$omp&firstprivate(wj)
!$omp do
do k=1,ysize(3)
  do i=1,ysize(1)
    call vrfftf(1,jtot,py(i,1:jtot,k),dj,1,wj)
  enddo
enddo
!$omp end parallel
!
call transpose_y_to_z(py,pz)
!
!$omp parallel default(none) &
!$omp&shared(pz,d,a,b,c,xyrt) &
!$omp&private(i,j,k,z,bb)
!$omp do
do j=1,jmax
  do i=1,imax
    z        = 1./( bbb(1) - coef*xyrt(i,j) )
    d(i,j,1) = ccc(1)*z
    pz(i,j,1) = pz(i,j,1)*z
  enddo
enddo
!$omp barrier
do k=2,kkmax-1
!$omp do
   do j=1,jmax
     do i=1,imax
       bb       = bbb(k) - coef*xyrt(i,j)
       z        = 1./ ( bb - aaa(k)*d(i,j,k-1) )
       d(i,j,k) = ccc(k)*z
       pz(i,j,k) = (pz(i,j,k) - aaa(k)*pz(i,j,k-1))*z
     enddo
  enddo
!$omp barrier
enddo
!$omp do
do j=1,jmax
  do i=1,imax
    bb       = bbb(kkmax) - coef*xyrt(i,j)
    z        = bb - aaa(kkmax)*d(i,j,kkmax-1)
    if(z.ne.0.) then
      pz(i,j,kkmax) = (pz(i,j,kkmax) - aaa(kkmax)*pz(i,j,kkmax-1))/z
    else
      pz(i,j,kkmax) =0. 
    endif
  enddo
enddo

do k=kkmax-1,1,-1
!$omp do
  do j=1,jmax
    do i=1,imax
      pz(i,j,k) = pz(i,j,k)-d(i,j,k)*pz(i,j,k+1)
    enddo
  enddo
!$omp barrier
enddo
!$omp end parallel

call transpose_z_to_y(pz,py)
!$omp parallel default(none) &
!$omp& shared(py,ysize) private(i,j,k,dj) &
!$omp&firstprivate(wj)
!$omp do
do k=1,ysize(3)
  do i=1,ysize(1)
    call vrfftb(1,jtot,py(i,1:jtot,k),dj,1,wj)
  enddo
enddo
!$omp end parallel
!
call transpose_y_to_x(py,px)
!
!$omp parallel default(none) &
!$omp& shared(px,xsize) private(i,j,k,di) &
!$omp&firstprivate(wi)
!$omp do
do k=1,xsize(3)
  do j=1,xsize(2)
    call vrfftb(1,itot,px(1:itot,j,k),di,1,wi)
  enddo
enddo
!$omp end parallel
!
call transpose_x_to_y(px,py)
call transpose_y_to_z(py,pz)

return
end subroutine solver2d_Star


!
end module mod_solver
