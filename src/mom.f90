module mod_mom
use mod_param
use decomp_2d
use decomp_2d_io
use mpi
use mod_common_mpi
use mod_bound
use mod_interface
use mod_common
implicit none
private
public momxad,momxp,momyad,momyp,momzad,momzp,PFM,PFM_couple_x,PFM_couple_y,PFM_couple_z,momdiffvof, &
        momad,PFM_couple, momdiff,LapVel
contains
!
subroutine momxad(advnU,dfunU,u,v,w)
implicit none
real, dimension(0:,0:,0:), intent(in) :: u,v,w
real, dimension(0:,0:,0:), intent(out) :: advnU,dfunU
integer ::ipp,jpp,kpp, im,ip,jm,jp,km,kp,i,j,k
real :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
real :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
real:: dudx,dudy,dudz,dVisdx,dVisdy,dVisdz
real:: visci,viscip,viscim,viscjp,viscjm,visckp,visckm
real:: temp
real::rhoi,visi
real::g_x


do k=1,kmax
  do j=1,jmax
    do i=1,imax
      ip = i + 1 
      jp = j + 1 
      kp = k + 1 
      im = i - 1 
      jm = j - 1 
      km = k - 1 
      uuip  =( 0.25 * (U(i,j,k)+U(ip,j,k))*( U(i,j,k)+U(ip,j,k) ))* rhol(ip,j,k)
      uuim  = (0.25 * (U(im,j,k)+U(i,j,k))*( U(im,j,k)+U(i,j,k) ))* rhol(i,j,k)
      uvjp  =( 0.25 * (U(i,j,k )+U(i,jp,k) )*( V(i,j,k)+V(ip,j,k) ))* &
               0.25*(rhol(i,j,k)+rhol(ip,j,k)+rhol(ip,jp,k)+rhol(i,jp,k))
      uvjm  =( 0.25 * (U(i,j,k )+U(i,jm,k) )*( V(i,jm,k)+V(ip,jm,k) ))* &
               0.25*(rhol(i,j,k)+rhol(ip,j,k)+rhol(ip,jm,k)+rhol(i,jm,k))
      uwkp  =( 0.25 * (U(i,j,k )+U(i,j,kp) )*( W(i,j,k)+W(ip,j,k) ))* &
               0.25*(rhol(i,j,k)+rhol(ip,j,k)+rhol(ip,j,kp)+rhol(i,j,kp))
      uwkm  =( 0.25 * (U(i,j,k )+U(i,j,km) )*( W(i,j,km)+W(ip,j,km) ))* &
               0.25*(rhol(i,j,k)+rhol(ip,j,k)+rhol(ip,j,km)+rhol(i,j,km))

      dudxp = (U(ip,j,k)-U(i,j,k))*dxi
      dudxm = (U(i,j,k)-U(im,j,k))*dxi
      dudyp = (U(i,jp,k)-U(i,j,k))*dyi
      dudym = (U(i,j,k)-U(i,jm,k))*dyi
      dudzp = (U(i,j,kp)-U(i,j,k))*dzi
      dudzm = (U(i,j,k)-U(i,j,km))*dzi

      rhoi=0.5*(rhol(ip,j,k)+rhol(i,j,k))
      visi=0.5*(visl(ip,j,k)+visl(i,j,k))


      dfunU(i,j,k) = (visi/rhoi)* &
                     ((dudxp-dudxm)*dxi+(dudyp-dudym)*dyi+(dudzp-dudzm)*dzi)+g_x !+ &
                    !(1./rhoj)*(dvdx*dVisdx+dvdy*dVisdy+dvdz*dVisdz)+g_y 

      advnU(i,j,k) = (1./rhoi)*(dxi*( -uuip + uuim ) + & 
                     dyi*( -uvjp + uvjm ) + & 
                     dzi*( -uwkp + uwkm ) ) 



    enddo
  enddo
enddo
!$omp end parallel
return
end subroutine momxad



subroutine momxp(dudt,p)
implicit none
real, dimension(0:,0:,0:), intent(in) :: p
real, dimension(0:,0:,0:), intent(out) :: dudt
integer :: i,j,k,ip
!
do k=1,kmax
  do j=1,jmax
    do i=1,imax
        ip = i + 1 
        dudt(i,j,k) = -dxi*( p(ip,j,k)-p(i,j,k) ) *2.0/(rhol(ip,j,k)+rhol(i,j,k))
    enddo
  enddo
enddo
!
return
end subroutine momxp

!
subroutine momyad(advnV,dfunV,u,v,w)
implicit none
real, dimension(0:,0:,0:), intent(in) :: u,v,w
real, dimension(0:,0:,0:), intent(out) :: advnV,dfunV
integer :: ipp,jpp,kpp,im,ip,jm,jp,km,kp,i,j,k
real :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
real :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
real :: leftbound,rightbound
real:: dvdx,dvdy,dvdz,dVisdx,dVisdy,dVisdz
real:: visci,viscip,viscim,viscjp,viscjm,visckp,visckm
real::temp
real::rhoj,visj
real:: uvi,vvi,vwi,drhod1,drhod2,drhodx,drhody,drhodz

do k=1,kmax
  do j=1,jmax
    do i=1,imax
      ip = i + 1
      jp = j + 1
      kp = k + 1
      im = i - 1
      jm = j - 1
      km = k - 1
        uvip  =( 0.25 * (U(i,j,k)+U(i,jp,k))*( V(i,j,k)+V(ip,j,k) ))* &
                 0.25* (rhol(i,j,k)+rhol(ip,j,k)+rhol(ip,jp,k)+rhol(i,jp,k))
        uvim  = (0.25 * (U(im,j,k)+U(im,jp,k))*( V(i,j,k)+V(im,j,k) ))* &
                 0.25* (rhol(i,j,k)+rhol(im,j,k)+rhol(im,jp,k)+rhol(i,jp,k))
        vvjp  =( 0.25 * (V(i,j,k )+V(i,jp,k) )*( V(i,j,k)+V(i,jp,k) ))*rhol(i,jp,k)
        vvjm  = 0.25 * (V(i,j,k )+V(i,jm,k) )*( V(i,j,k)+V(i,jm,k) )*rhol(i,j,k)
        wvkp  = (0.25 * (W(i,j,k )+W(i,jp,k) )*( V(i,j,kp)+V(i,j,k)) )* &
                 0.25* (rhol(i,j,k)+rhol(i,j,kp)+rhol(i,jp,kp)+rhol(i,jp,k))
        wvkm  = 0.25 * (W(i,j,km)+W(i,jp,km))*( V(i,j,km)+V(i,j,k) )* &
                0.25* (rhol(i,j,k)+rhol(i,j,km)+rhol(i,jp,km)+rhol(i,jp,k))
      dvdxp = (V(ip,j,k)-V(i,j,k))*dxi
      dvdxm = (V(i,j,k)-V(im,j,k))*dxi
      dvdyp = (V(i,jp,k)-V(i,j,k))*dyi
      dvdym = (V(i,j,k)-V(i,jm,k))*dyi
      dvdzp = (V(i,j,kp)-V(i,j,k))*dzi
      dvdzm = (V(i,j,k)-V(i,j,km))*dzi

      rhoj=0.5*(rhol(i,jp,k)+rhol(i,j,k))
      visj=0.5*(visl(i,jp,k)+visl(i,j,k))


      dfunV(i,j,k) = (visj/rhoj)* &
                     ((dvdxp-dvdxm)*dxi+(dvdyp-dvdym)*dyi+(dvdzp-dvdzm)*dzi)+g_y !+ &
                    !(1./rhoj)*(dvdx*dVisdx+dvdy*dVisdy+dvdz*dVisdz)+g_y 

      advnV(i,j,k) = (1./rhoj)*(dxi*( -uvip + uvim ) + &
                     dyi*( -vvjp + vvjm ) + &
                     dzi*( -wvkp + wvkm ) )



    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momyad
!
subroutine momyp(dvdt,p)
implicit none
real, dimension(0:,0:,0:), intent(in) :: p
real, dimension(0:,0:,0:), intent(out) :: dvdt
integer :: i,j,k,jp
!
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      jp = j + 1
        dvdt(i,j,k) = - dyi*( p(i,jp,k)-p(i,j,k) ) * 2.0/(rhol(i,jp,k)+rhol(i,j,k))
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momyp
!
subroutine momzad(advnW,dfunW,u,v,w)
implicit none
real, dimension(0:,0:,0:), intent(in) :: u,v,w
real, dimension(0:,0:,0:), intent(out) :: advnW,dfunW
integer ::ipp,jpp,kpp, im,ip,jm,jp,km,kp,i,j,k
real :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
real :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
real :: leftbound,rightbound
real:: dwdx,dwdy,dwdz,dVisdx,dVisdy,dVisdz
real:: visci,viscip,viscim,viscjp,viscjm,visckp,visckm
real::visk,rhok
real:: uwi,vwi,wwi,drhod1,drhod2,drhodx,drhody,drhodz

do k=1,kmax
  do j=1,jmax
    do i=1,imax
      ip = i + 1
      jp = j + 1
      kp = k + 1
      im = i - 1
      jm = j - 1
      km = k - 1

!if (2.gt.3) then

      uwip  = (0.25 * ( W(i,j,k)+W(ip,j,k))*(U(i ,j,k)+U(i ,j,kp) ))* &
              0.25* (rhol(i,j,k)+rhol(ip,j,k)+rhol(ip,j,kp)+rhol(i,j,kp))


      uwim  = (0.25 * ( W(i,j,k)+W(im,j,k))*(U(im,j,k)+U(im,j,kp) ))* &
              0.25* (rhol(i,j,k)+rhol(i,j,kp)+rhol(im,j,k)+rhol(im,j,kp))



      vwjp  = (0.25 * ( W(i,j,k)+W(i,jp,k))*(V(i ,j,k)+V(i,j ,kp) ))*&
            0.25* (rhol(i,j,k)+rhol(i,jp,k)+rhol(i,jp,kp)+rhol(i,j,kp))


      vwjm  = (0.25 * ( W(i,j,k)+W(i,jm,k))*(V(i ,jm,k)+V(i,jm ,kp) ))*&
           0.25* (rhol(i,j,k)+rhol(i,j,kp)+rhol(i,jm,kp)+rhol(i,jm,k))


      wwkp  = (0.25 * ( W(i,j,k)+W(i,j,kp))*(W(i ,j,k)+W(i,j,kp ) )) * rhol(i,j,kp)



      wwkm  = (0.25 * ( W(i,j,k)+W(i,j,km))*(W(i ,j,k)+W(i,j,km ) ))*rhol(i,j,k)




      dwdxp = (W(ip,j,k)-W(i,j,k))*dxi
      dwdxm = (W(i,j,k)-W(im,j,k))*dxi
      dwdyp = (W(i,jp,k)-W(i,j,k))*dyi
      dwdym = (W(i,j,k)-W(i,jm,k))*dyi
      dwdzp = (W(i,j,kp)-W(i,j,k))*dzi
      dwdzm = (W(i,j,k)-W(i,j,km))*dzi



      rhok= 0.5* (rhol(i,j,kp)+rhol(i,j,k))
      visk= 0.5* (visl(i,j,kp)+visl(i,j,k))
      
      advnW(i,j,k) =(1./rhok)*( dxi*( -uwip + uwim ) + &
                     dyi*( -vwjp + vwjm ) + &
                     dzi*( -wwkp + wwkm ))
                       


      dfunW(i,j,k) = (visk/rhok)* &
                     ((dwdxp-dwdxm)*dxi+(dwdyp-dwdym)*dyi+(dwdzp-dwdzm)*dzi)+g_z !+ &
                    ! (1./rhok)*(dwdx*dVisdx+dwdy*dVisdy+dwdz*dVisdz)+g_z 

    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momzad
!
subroutine momzp(dwdt,p)
implicit none
real, dimension(0:,0:,0:), intent(in) :: p
real, dimension(0:,0:,0:), intent(out) :: dwdt
integer :: i,j,k,kp
!
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      kp = k + 1
        dwdt(i,j,k) = - dzi*( p(i,j,kp)-p(i,j,k) ) * 2.0/(rhol(i,j,kp)+rhol(i,j,k))
    enddo
  enddo
enddo
!
return
end subroutine momzp

!****************************** Phase Field Model *****************************************************
subroutine  PFM(dPFM,PFM_phi)
implicit none
real, dimension(-5:i1+5,-5:j1+5,-5:k1+5), intent(in) :: PFM_phi
real, dimension(0:,0:,0:), intent(out) :: dPFM
real:: RHS_all
integer:: i,j,k
!!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(uvip,uvim,vvjp,vvjm,wvkp,wvkm,dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm)
!!$omp do


do k=1,kmax
  do j=1,jmax
    do i=1,imax
dPFM(i,j,k)=0.0
call PFM_RHS(RHS_all,PFM_phi,i,j,k)
dPFM(i,j,k)=RHS_all
    enddo
  enddo
enddo
!!$omp end parallel
!
return
end subroutine PFM
!
subroutine  PFM_couple_x(dPFM_couple_x)
implicit none
real, dimension(0:,0:,0:), intent(out) ::dPFM_couple_x
real::Fx_SurfTen,Fy_SurfTen,Fz_SurfTen
integer::i,j,k
do k=1,kmax
  do j=1,jmax
    do i=1,imax
       dPFM_couple_x(i,j,k)=0.0
       call PFM_SurfaceTension (i,j,k,Fx_SurfTen,Fy_SurfTen,Fz_SurfTen)
       dPFM_couple_x(i,j,k)= Fx_SurfTen
    enddo
  enddo
enddo
return
end subroutine PFM_couple_x


subroutine  PFM_couple_y(dPFM_couple_y)
implicit none
real, dimension(0:,0:,0:), intent(out) ::dPFM_couple_y
integer:: i,j,k
real::Fx_SurfTen,Fy_SurfTen,Fz_SurfTen

!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(uvip,uvim,vvjp,vvjm,wvkp,wvkm,dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm)
!$omp do


do k=1,kmax
  do j=1,jmax
    do i=1,imax
      dPFM_couple_y(i,j,k)=0.0
      call PFM_SurfaceTension (i,j,k,Fx_SurfTen,Fy_SurfTen,Fz_SurfTen)
      dPFM_couple_y(i,j,k) = Fy_SurfTen
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine PFM_couple_y

subroutine  PFM_couple_z(dPFM_couple_z)
implicit none
real, dimension(0:,0:,0:), intent(out) ::dPFM_couple_z
integer:: i,j,k
real::Fx_SurfTen,Fy_SurfTen,Fz_SurfTen
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(uvip,uvim,vvjp,vvjm,wvkp,wvkm,dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm)
!$omp do
do k=1,kmax
  do j=1,jmax
    do i=1,imax
     call PFM_SurfaceTension (i,j,k,Fx_SurfTen,Fy_SurfTen,Fz_SurfTen)
     dPFM_couple_z(i,j,k) = Fz_SurfTen
    enddo
  enddo
enddo
!$omp end parallel
return
end subroutine PFM_couple_z


subroutine momdiffvof(duDiffVof,dvDiffVof,dwDiffVof,u,v,w,C)
implicit none
real, dimension(0:i1,0:j1,0:k1), intent(in) :: u,v,w
real, dimension(-5:i1+5,-5:j1+5,-5:k1+5), intent(in) :: C
real, dimension(0:i1,0:j1,0:k1), intent(out) :: duDiffVof, dvDiffVof, dwDiffVof
integer ::  i,j,k,im,ip,jm,jp,km,kp
real :: VGRux,VGRuy,VGRuz,VGRvx,VGRvy,VGRvz,VGRwx,VGRwy,VGRwz
real :: GRx11L,GRy11L,GRz11L
real:: phi,rhoijk
do k=1,kmax
  do j=1,jmax
    do i=1,imax
! diffusion

      ip = i + 1
      jp = j + 1
      kp = k + 1
      im = i - 1
      jm = j - 1
      km = k - 1

      VGRux = 0.5*(U(ip,j,k) + U(i,j,k))*dxi - 0.5*(U(i,j,k) + U(im,j,k))*dxi
      VGRuy = 0.5*(U(i,jp,k) + U(i,j,k))*dyi - 0.5*(U(i,j,k) + U(i,jm,k))*dyi
      VGRuz = 0.5*(U(i,j,kp) + U(i,j,k))*dzi - 0.5*(U(i,j,k) + U(i,j,km))*dzi
      VGRvx = 0.5*(V(ip,j,k) + V(ip,jm,k))*dxi - 0.5*(V(i,j,k) + V(i,jm,k))*dxi
      VGRwx = 0.5*(W(ip,j,k) + W(ip,j,km))*dxi - 0.5*(W(i,j,k) + W(i,j,km))*dxi

      GRx11L = (C(ip,j,k) - C(i,j,k))*dxi
      GRy11L = 0.25*(C(ip,jp,k) + C(ip,j,k) + C(i,jp,k) + C(i,j,k))*dyi - &
               0.25*(C(ip,jm,k) + C(ip,j,k) + C(i,jm,k) + C(i,j,k))*dyi
      GRz11L = 0.25*(C(ip,j,kp) + C(i,j,kp) + C(ip,j,k) + C(i,j,k))*dzi - &
               0.25*(C(ip,j,km) + C(i,j,km) + C(ip,j,k) + C(i,j,k))*dzi

      phi = 0.5*(C(i,j,k) + C(ip,j,k))
      rhoijk= 0.5*(phi+1)*rho_2-0.5*(phi-1)*rho_1
      duDiffVof(i,j,k) = 0.5*(1./rhoijk)*(vis_2-vis_1)*(GRx11L*(VGRux+VGRux)+ &
                                      GRy11L*(VGRuy+VGRvx)+ &
                                      GRz11L*(VGRuz+VGRwx))

      VGRvx = 0.5*(V(ip,j,k) + V(i,j,k))*dxi - 0.5*(V(i,j,k) + V(im,j,k))*dxi
      VGRvy = 0.5*(V(i,jp,k) + V(i,j,k))*dyi - 0.5*(V(i,j,k) + V(i,jm,k))*dyi
      VGRvz = 0.5*(V(i,j,kp) + V(i,j,k))*dzi - 0.5*(V(i,j,k) + V(i,j,km))*dzi
      VGRuy = 0.5*(U(i,jp,k) + U(im,jp,k))*dyi - 0.5*(U(i,j,k) + U(im,j,k))*dyi
      VGRwy = 0.5*(W(i,jp,k) + W(i,jp,km))*dyi - 0.5*(W(i,j,k) + W(i,j,km))*dyi

      GRx11L = 0.25*(C(ip,jp,k) + C(ip,j,k) + C(i,jp,k) + C(i,j,k))*dxi - &
               0.25*(C(im,jp,k) + C(im,j,k) + C(i,jp,k) + C(i,j,k))*dxi
      GRy11L = (C(i,jp,k) - C(i,j,k))*dyi
      GRz11L = 0.25*(C(i,jp,kp) + C(i,j,kp) + C(i,jp,k) + C(i,j,k))*dzi - &
               0.25*(C(i,jp,km) + C(i,j,km) + C(i,jp,k) + C(i,j,k))*dzi

      phi = 0.5*(C(i,j,k) + C(i,jp,k))
      rhoijk= 0.5*(phi+1)*rho_2-0.5*(phi-1)*rho_1
      dvDiffVof(i,j,k) = 0.5*(1./rhoijk)*(vis_2-vis_1)*(GRx11L*(VGRvx+VGRuy)+ &
                                      GRy11L*(VGRvy+VGRvy)+ &
                                      GRz11L*(VGRvz+VGRwy))

      VGRwx = 0.5*(W(ip,j,k) + W(i,j,k))*dxi - 0.5*(W(i,j,k) + W(im,j,k))*dxi
      VGRwy = 0.5*(W(i,jp,k) + W(i,j,k))*dyi - 0.5*(W(i,j,k) + W(i,jm,k))*dyi
      VGRwz = 0.5*(W(i,j,kp) + W(i,j,k))*dzi - 0.5*(W(i,j,k) + W(i,j,km))*dzi
      VGRuz = 0.5*(U(i,j,kp) + U(im,j,kp))*dzi - 0.5*(U(i,j,k) + U(im,j,k))*dzi
      VGRvz = 0.5*(V(i,j,kp) + V(i,jm,kp))*dzi - 0.5*(V(i,j,k) + V(i,jm,k))*dzi

      GRx11L = 0.25*(C(ip,j,kp) + C(i,j,kp) + C(ip,j,k) + C(i,j,k))*dxi - &
               0.25*(C(im,j,kp) + C(i,j,kp) + C(im,j,k) + C(i,j,k))*dxi
      GRy11L = 0.25*(C(i,jp,kp) + C(i,j,kp) + C(i,jp,k) + C(i,j,k))*dyi - &
               0.25*(C(i,jm,kp) + C(i,j,kp) + C(i,jm,k) + C(i,j,k))*dyi
      GRz11L = (C(i,j,kp) - C(i,j,k))*dzi

      phi = 0.5*(C(i,j,k) + C(i,j,kp))
      rhoijk= 0.5*(phi+1)*rho_2-0.5*(phi-1)*rho_1
      dwDiffVof(i,j,k) = 0.5*(1./rhoijk)*(vis_2-vis_1)*(GRx11L*(VGRwx+VGRuz)+ &
                                      GRy11L*(VGRwy+VGRvz)+ &
                                      GRz11L*(VGRwz+VGRwz))

    enddo
  enddo
enddo
!
return
end subroutine momdiffvof

subroutine momad( advnU,advnV,advnW,u,v,w)
implicit none
real, dimension(0:i1,0:j1,0:k1), intent(in) :: u,v,w
real, dimension(0:i1,0:j1,0:k1), intent(out) :: advnU, advnV, advnW
real :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
real :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
real :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
real :: leftbound,rightbound
integer::i,j,k,ip,jp,kp,im,jm,km
 
do k=1,kmax
  do j=1,jmax
    do i=1,imax
         ip = i+1 
         im = i-1 
         jp = j+1 
         jm = j-1 
         kp = k+1 
         km = k-1 

      ! U-component
         uuip  = 0.25 * ( U(ip,j,k)+U(i,j,k) )*( U(ip,j,k)+U(i,j,k)  )
         uuim  = 0.25 * ( U(im,j,k)+U(i,j,k) )*( U(im,j,k)+U(i,j,k)  )
         uvjp  = 0.25 * ( U(i,jp,k)+U(i,j,k) )*( V(ip,j,k)+V(i,j,k)  )
         uvjm  = 0.25 * ( U(i,jm,k)+U(i,j,k) )*( V(ip,jm,k)+V(i,jm,k))
         uwkp  = 0.25 * ( U(i,j,kp)+U(i,j,k) )*( W(ip,j,k) +W(i,j,k) )
         uwkm  = 0.25 * ( U(i,j,km)+U(i,j,k) )*( W(ip,j,km)+W(i,j,km))
         advnU(i,j,k) = dxi*( -uuip + uuim ) + & 
                     dyi*( -uvjp + uvjm ) + & 
                     dzi*( -uwkp + uwkm ) + g_x 

      ! V-component
        uvip  = 0.25 * (U(i,j,k) +U(i,jp,k)) *( V(i,j,k) +V(ip,j,k) )
        uvim  = 0.25 * (U(im,j,k)+U(im,jp,k))*( V(i,j,k) +V(im,j,k) )
        vvjp  = 0.25 * (V(i,j,k )+V(i,jp,k) )*( V(i,j,k) +V(i,jp,k) )
        vvjm  = 0.25 * (V(i,j,k )+V(i,jm,k) )*( V(i,j,k) +V(i,jm,k) )
        wvkp  = 0.25 * (W(i,j,k )+W(i,jp,k) )*( V(i,j,kp)+V(i,j,k)  )   
        wvkm  = 0.25 * (W(i,j,km)+W(i,jp,km))*( V(i,j,km)+V(i,j,k)  )
        advnV(i,j,k) =  dxi*( -uvip + uvim ) + & 
                        dyi*( -vvjp + vvjm ) + & 
                        dzi*( -wvkp + wvkm ) + g_y 

      ! W-Component
        uwip  = 0.25 * ( W(i,j,k)+W(ip,j,k))*(U(i ,j,k) +U(i ,j,kp) )
        uwim  = 0.25 * ( W(i,j,k)+W(im,j,k))*(U(im,j,k) +U(im,j,kp) )
        vwjp  = 0.25 * ( W(i,j,k)+W(i,jp,k))*(V(i ,j,k) +V(i,j ,kp) )
        vwjm  = 0.25 * ( W(i,j,k)+W(i,jm,k))*(V(i ,jm,k)+V(i,jm ,kp))
        wwkp  = 0.25 * ( W(i,j,k)+W(i,j,kp))*(W(i ,j,k) +W(i,j,kp)  )   
        wwkm  = 0.25 * ( W(i,j,k)+W(i,j,km))*(W(i ,j,k) +W(i,j,km)  )
        advnW(i,j,k) = dxi*( -uwip + uwim ) + & 
                       dyi*( -vwjp + vwjm ) + & 
                       dzi*( -wwkp + wwkm ) +g_z

    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momad
subroutine  PFM_couple(u,v,w,dPFM_couple_x,dPFM_couple_y,dPFM_couple_z)
implicit none
real, dimension(0:i1,0:j1,0:k1), intent(out) ::dPFM_couple_x,dPFM_couple_y,dPFM_couple_z
real, dimension(0:i1,0:j1,0:k1), intent(in)   :: u,v,w
real:: phiijk,phiip,phiim,phijp,phijm,phikp,phikm,grad_phi
integer::i,j,k,ip,im,jp,jm,kp,km
real::rhoijk ,Psi_ijk
real:: vel_ip,vel_im,vel_jp,vel_jm,vel_kp,vel_km
real::chem_ip,chem_im,chem_jp,chem_jm,chem_kp,chem_km,high_dens
do k=1,kmax
  do j=1,jmax
    do i=1,imax
     ip = i+1 
     im = i-1 
     jp = j+1 
     jm = j-1 
     kp = k+1 
     km = k-1 
! ------ U-faces------------------------
     dPFM_couple_x(i,j,k)=0.0
     phiijk        = 0.5*(PFM_phi(i,j,k)+PFM_phi(ip,j,k))
     rhoijk       = 0.5*(phiijk+1)*rho_2-0.5*(phiijk-1)*rho_1
     phiip         = PFM_phi(ip,j,k)
     phiim         = PFM_phi(i,j,k)
     phijp         = 0.25*(PFM_phi(i,j,k)+PFM_phi(ip,j,k)+PFM_phi(ip,jp,k)+PFM_phi(i,jp,k))
     phijm         = 0.25*(PFM_phi(i,j,k)+PFM_phi(ip,j,k)+PFM_phi(ip,jm,k)+PFM_phi(i,jm,k))
     phikp         = 0.25*(PFM_phi(i,j,k)+PFM_phi(ip,j,k)+PFM_phi(ip,j,kp)+PFM_phi(i,j,kp))
     phikm         = 0.25*(PFM_phi(i,j,k)+PFM_phi(ip,j,k)+PFM_phi(ip,j,km)+PFM_phi(i,j,km))
     grad_phi      = (phiip-phiim)*dxi
     vel_ip               = 0.5*(u(i,j,k)+u(i+1,j,k))
     vel_im               = 0.5*(u(i,j,k)+u(i-1,j,k))
     vel_jp               = 0.5*(u(i,j,k)+u(i,j+1,k))
     vel_jm               = 0.5*(u(i,j,k)+u(i,j-1,k))
     vel_kp               = 0.5*(u(i,j,k)+u(i,j,k+1))
     vel_km               = 0.5*(u(i,j,k)+u(i,j,k-1))
     chem_ip              = chem_pot(i+1,j,k)
     chem_im              = chem_pot(i,j,k)
     chem_jp              = 0.25*(chem_pot(i,j,k)+chem_pot(i+1,j,k)+chem_pot(i+1,j+1,k)+chem_pot(i,j+1,k))
     chem_jm              = 0.25*(chem_pot(i,j,k)+chem_pot(i+1,j,k)+chem_pot(i+1,j-1,k)+chem_pot(i,j-1,k))
     chem_kp              = 0.25*(chem_pot(i,j,k)+chem_pot(i+1,j,k)+chem_pot(i+1,j,k+1)+chem_pot(i,j,k+1))
     chem_km              = 0.25*(chem_pot(i,j,k)+chem_pot(i+1,j,k)+chem_pot(i+1,j,k-1)+chem_pot(i,j,k-1))
     high_dens            = 0.5*(rho_2-rho_1)*mobility*( (chem_ip-chem_im)*dxi*(vel_ip-vel_im)*dxi+ &
                                                         (chem_jp-chem_jm)*dyi*(vel_jp-vel_jm)*dyi+ &
                                                         (chem_kp-chem_km)*dzi*(vel_kp-vel_km)*dzi)
     dPFM_couple_x(i,j,k) = 0.5*(chem_pot(i,j,k)+chem_pot(ip,j,k))*grad_phi/rhoijk !+ high_dens/rhoijk
! ------ V-faces------------------------
     dPFM_couple_y(i,j,k)=0.0
     phiijk               = 0.5*(PFM_phi(i,j,k)+PFM_phi(i,jp,k))
     rhoijk               = 0.5*(phiijk+1)*rho_2-0.5*(phiijk-1)*rho_1
     phijp                = PFM_phi(i,jp,k)
     phijm                = PFM_phi(i,j,k)
     grad_phi             = (phijp-phijm)*dyi
     vel_ip               = 0.5*(v(i,j,k)+v(i+1,j,k))
     vel_im               = 0.5*(v(i,j,k)+v(i-1,j,k))
     vel_jp               = 0.5*(v(i,j,k)+v(i,j+1,k))
     vel_jm               = 0.5*(v(i,j,k)+v(i,j-1,k))
     vel_kp               = 0.5*(v(i,j,k)+v(i,j,k+1))
     vel_km               = 0.5*(v(i,j,k)+v(i,j,k-1))
     chem_ip              = 0.25*(chem_pot(i,j,k)+chem_pot(i,j+1,k)+chem_pot(i+1,j+1,k)+chem_pot(i+1,j,k))
     chem_im              = 0.25*(chem_pot(i,j,k)+chem_pot(i,j+1,k)+chem_pot(i-1,j+1,k)+chem_pot(i-1,j,k))
     chem_jp              = chem_pot(i,j+1,k)
     chem_jm              = chem_pot(i,j,k)
     chem_kp              = 0.25*(chem_pot(i,j,k)+chem_pot(i,j+1,k)+chem_pot(i,j+1,k+1)+chem_pot(i,j,k+1))
     chem_km              = 0.25*(chem_pot(i,j,k)+chem_pot(i,j+1,k)+chem_pot(i,j+1,k-1)+chem_pot(i,j,k-1))
     high_dens            = 0.5*(rho_2-rho_1)*mobility*( (chem_ip-chem_im)*dxi*(vel_ip-vel_im)*dxi+ &
                                                         (chem_jp-chem_jm)*dyi*(vel_jp-vel_jm)*dyi+ &
                                                         (chem_kp-chem_km)*dzi*(vel_kp-vel_km)*dzi)

     dPFM_couple_y(i,j,k) = 0.5*(chem_pot(i,j,k)+chem_pot(i,jp,k))*grad_phi/rhoijk !+ high_dens/rhoijk 
 ! ------ w-faces------------------------
     dPFM_couple_z(i,j,k)=0.0
     phiijk               = 0.5*(PFM_phi(i,j,k)+PFM_phi(i,j,kp))
     rhoijk               = 0.5*(phiijk+1)*rho_2-0.5*(phiijk-1)*rho_1
     phikp                = PFM_phi(i,j,kp)
     phikm                = PFM_phi(i,j,k)
     grad_phi             = (phikp-phikm)*dzi
     vel_jp               = 0.5*(w(i,j,k)+w(i+1,j,k))
     vel_jm               = 0.5*(w(i,j,k)+w(i-1,j,k))
     vel_jp               = 0.5*(w(i,j,k)+w(i,j+1,k))
     vel_jm               = 0.5*(w(i,j,k)+w(i,j-1,k))
     vel_kp               = 0.5*(w(i,j,k)+w(i,j,k+1))
     vel_km               = 0.5*(w(i,j,k)+w(i,j,k-1))
     chem_ip              = 0.25*(chem_pot(i,j,k)+chem_pot(i+1,j,k)+chem_pot(i+1,j,k+1)+chem_pot(i,j,k+1))
     chem_im              = 0.25*(chem_pot(i,j,k)+chem_pot(i-1,j,k)+chem_pot(i-1,j,k+1)+chem_pot(i,j,k+1))
     chem_jp              = 0.25*(chem_pot(i,j,k)+chem_pot(i,j+1,k)+chem_pot(i,j+1,k+1)+chem_pot(i,j,k+1))
     chem_jm              = 0.25*(chem_pot(i,j,k)+chem_pot(i,j-1,k)+chem_pot(i,j-1,k+1)+chem_pot(i,j,k+1))
     chem_kp              = chem_pot(i,j,k+1)
     chem_km              = chem_pot(i,j,k)
     high_dens            = 0.5*(rho_2-rho_1)*mobility*( (chem_ip-chem_im)*dxi*(vel_ip-vel_im)*dxi+ &
                                                         (chem_jp-chem_jm)*dyi*(vel_jp-vel_jm)*dyi+ &
                                                         (chem_kp-chem_km)*dzi*(vel_kp-vel_km)*dzi)


      dPFM_couple_z(i,j,k) = 0.5*(chem_pot(i,j,k)+chem_pot(i,j,kp))*grad_phi/rhoijk !+ high_dens/rhoijk
! ------ w-faces------------------------
    enddo
  enddo
enddo
return
end subroutine PFM_couple

subroutine LapVel(lapU,lapV,lapW,u,v,w)
implicit none
real, dimension(0:i1,0:j1,0:k1), intent(in) :: u,v,w
real, dimension(0:i1,0:j1,0:k1), intent(out) :: lapU,lapV,lapW
integer :: im,ip,jm,jp,km,kp,i,j,k
real :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
real :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
real :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm


do k=1,kmax
  do j=1,jmax
    do i=1,imax
      ip = i + 1
      jp = j + 1
      kp = k + 1
      im = i - 1
      jm = j - 1
      km = k - 1
      dudxp = (U(ip,j,k)-U(i,j,k))*dxi
      dudxm = (U(i,j,k)-U(im,j,k))*dxi
      dudyp = (U(i,jp,k)-U(i,j,k))*dyi
      dudym = (U(i,j,k)-U(i,jm,k))*dyi
      dudzp = (U(i,j,kp)-U(i,j,k))*dzi
      dudzm = (U(i,j,k)-U(i,j,km))*dzi
      lapU(i,j,k) = (dudxp-dudxm)*dxi + &
                     (dudyp-dudym)*dyi + &
                     (dudzp-dudzm)*dzi


      dvdxp = (V(ip,j,k)-V(i,j,k))*dxi
      dvdxm = (V(i,j,k)-V(im,j,k))*dxi
      dvdyp = (V(i,jp,k)-V(i,j,k))*dyi
      dvdym = (V(i,j,k)-V(i,jm,k))*dyi
      dvdzp = (V(i,j,kp)-V(i,j,k))*dzi
      dvdzm = (V(i,j,k)-V(i,j,km))*dzi
      lapV(i,j,k) = (dvdxp-dvdxm)*dxi + &
                     (dvdyp-dvdym)*dyi + &
                     (dvdzp-dvdzm)*dzi



      dwdxp = (W(ip,j,k)-W(i,j,k))*dxi
      dwdxm = (W(i,j,k)-W(im,j,k))*dxi
      dwdyp = (W(i,jp,k)-W(i,j,k))*dyi
      dwdym = (W(i,j,k)-W(i,jm,k))*dyi
      dwdzp = (W(i,j,kp)-W(i,j,k))*dzi
      dwdzm = (W(i,j,k)-W(i,j,km))*dzi
      lapW(i,j,k) = (dwdxp-dwdxm)*dxi + &
                     (dwdyp-dwdym)*dyi + &
                     (dwdzp-dwdzm)*dzi
    enddo
  enddo
enddo


end subroutine LapVel

subroutine momdiff(diffU,diffV,diffW,U,V,W)
implicit none
real, dimension(0:i1,0:j1,0:k1), intent(in) :: u,v,w
real, dimension(0:i1,0:j1,0:k1), intent(out) :: diffU,diffV,diffW
integer :: im,ip,jm,jp,km,kp,i,j,k
real:: visc_temp,duj_dxi,dui_dxj,Dx_ip,Dx_im,Dy_jp,Dy_jm
real::dDxdx,dDydy,Dz_kp,Dz_km,dDzdz,rhoijk
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      ip = i + 1
      jp = j + 1
      kp = k + 1
      im = i - 1
      jm = j - 1
      km = k - 1

   !*************** U-Component ******************************************
    ! dDx/dx
    visc_temp = visl(ip,j,k)
    duj_dxi   = (U(ip,j,k)-U(i,j,k))*dxi
    dui_dxj   = (U(ip,j,k)-U(i,j,k))*dxi
    Dx_ip     = visc_temp*(dui_dxj+duj_dxi)
    visc_temp = visl(i,j,k)
    duj_dxi   = (U(i,j,k)-U(im,j,k))*dxi
    dui_dxj   = (U(i,j,k)-U(im,j,k))*dxi
    Dx_im     =  visc_temp*(dui_dxj+duj_dxi)
    dDxdx     = (Dx_ip-Dx_im)*dxi
    ! dDy/dy
    visc_temp = 0.25*(visl(i,j,k)+visl(ip,j,k)+visl(ip,jp,k)+visl(i,jp,k))
    duj_dxi   = (U(i,jp,k)-U(i,j,k))*dyi
    dui_dxj   = (V(ip,j,k)-V(i,j,k))*dxi
    Dy_jp     =  visc_temp*(dui_dxj+duj_dxi)
    visc_temp =  0.25*(visl(i,j,k)+visl(ip,j,k)+visl(ip,jm,k)+visl(i,jm,k))
    duj_dxi   = (V(ip,jm,k)-V(i,jm,k))*dxi
    dui_dxj   = (U(i,j,k)-U(i,jm,k))*dyi
    Dy_jm     =  visc_temp*(dui_dxj+duj_dxi)
    dDydy     = (Dy_jp-Dy_jm)*dyi
    ! dDz/dz
    visc_temp = 0.25*(visl(i,j,k)+visl(ip,j,k)+visl(ip,j,kp)+visl(i,j,kp))
    duj_dxi   = (U(i,j,kp)-U(i,j,k))*dzi
    dui_dxj   = (W(ip,j,k)-W(i,j,k))*dxi
    Dz_kp     =  visc_temp*(dui_dxj+duj_dxi)
    visc_temp =  0.25*(visl(i,j,k)+visl(ip,j,k)+visl(ip,j,km)+visl(i,j,km))
    duj_dxi   = (W(ip,j,km)-W(i,j,km))*dxi
    dui_dxj   = (U(i,j,k)-U(i,j,km))*dzi
    Dz_km     =  visc_temp*(dui_dxj+duj_dxi)
    dDzdz     = (Dz_kp-Dz_km)*dzi
    rhoijk    = 0.5*(rhol(ip,j,k)+rhol(i,j,k))
    diffU(i,j,k) = (1./rhoijk)*(dDxdx+dDydy+dDzdz)
   !*************** V-Component ******************************************
    ! dDx/dx
    visc_temp = 0.25*(visl(i,j,k)+visl(i,jp,k)+visl(ip,jp,k)+visl(ip,j,k))
    dui_dxj   = (U(i,jp,k)-U(i,j,k))*dyi
    duj_dxi   = (V(ip,j,k)-V(i,j,k))*dxi
    Dx_ip      = visc_temp*(dui_dxj+duj_dxi)
    visc_temp = 0.25*(visl(i,j,k)+visl(i,jp,k)+visl(im,jp,k)+visl(im,j,k))
    dui_dxj   = (U(im,jp,k)-U(im,j,k))*dyi
    duj_dxi   = (V(i,j,k)-V(im,j,k))*dxi
    Dx_im      = visc_temp*(dui_dxj+duj_dxi)
    dDxdx     = (Dx_ip-Dx_im)*dxi
    ! dDy/dy
    visc_temp =  visl(i,jp,k)
    duj_dxi   = (V(i,jp,k)-V(i,j,k))*dyi
    dui_dxj   = (V(i,jp,k)-V(i,j,k))*dyi
    Dy_jp     =  visc_temp*(dui_dxj+duj_dxi)
    visc_temp =  visl(i,j,k)
    duj_dxi   = (V(i,j,k)-V(i,jm,k))*dyi
    dui_dxj   = (V(i,j,k)-V(i,jm,k))*dyi
    Dy_jm     =  visc_temp*(dui_dxj+duj_dxi)
    dDydy     = (Dy_jp-Dy_jm)*dyi
    ! dDz/dz
    visc_temp = 0.25*(visl(i,j,k)+visl(i,jp,k)+visl(i,jp,kp)+visl(i,j,kp))
    dui_dxj   = (W(i,jp,k)-W(i,j,k))*dyi
    duj_dxi   = (V(i,j,kp)-V(i,j,k))*dzi
    Dz_kp      = visc_temp*(dui_dxj+duj_dxi)
    visc_temp = 0.25*(visl(i,j,k)+visl(i,jp,k)+visl(i,jp,km)+visl(i,j,km))
    dui_dxj   = (W(i,jp,km)-W(i,j,km))*dyi
    duj_dxi   = (V(i,j,k)-V(i,j,km))*dzi
    Dz_km      = visc_temp*(dui_dxj+duj_dxi)
    dDzdz     = (Dz_kp-Dz_km)*dzi
    rhoijk    = 0.5*(rhol(i,j,k)+rhol(i,jp,k))
    diffV(i,j,k) = (1./rhoijk)*(dDxdx+dDydy+dDzdz)
   !*************** W-Component ******************************************
    ! dDx/dx
    visc_temp = 0.25*(visl(i,j,k)+visl(ip,j,k)+visl(ip,j,kp)+visl(i,j,kp))
    duj_dxi   = (U(i,j,kp)-U(i,j,k))*dzi
    dui_dxj   = (W(ip,j,k)-W(i,j,k))*dxi
    Dx_ip     =  visc_temp*(dui_dxj+duj_dxi)
    visc_temp =  0.25*(visl(i,j,k)+visl(im,j,k)+visl(im,j,kp)+visl(i,j,kp))
    duj_dxi   = (U(im,j,kp)-U(im,j,k))*dzi
    dui_dxj   = (W(i,j,k)-W(im,j,k))*dxi
    Dx_im     =  visc_temp*(dui_dxj+duj_dxi)
    dDxdx     = (Dx_ip-Dx_im)*dxi
    ! dDy/dy
    visc_temp = 0.25*(visl(i,j,k)+visl(i,jp,k)+visl(i,jp,kp)+visl(i,j,kp))
    duj_dxi   = (V(i,j,kp)-V(i,j,k))*dzi
    dui_dxj   = (W(i,jp,k)-W(i,j,k))*dyi
    Dy_jp     =  visc_temp*(dui_dxj+duj_dxi)
    visc_temp =  0.25*(visl(i,j,k)+visl(i,jm,k)+visl(i,jm,kp)+visl(i,j,kp))
    duj_dxi   = (V(i,jm,kp)-V(i,jm,k))*dzi
    dui_dxj   = (W(i,j,k)-W(i,jm,k))*dyi
    Dy_jm     =  visc_temp*(dui_dxj+duj_dxi)
    dDydy     = (Dy_jp-Dy_jm)*dyi
    ! dDz/dz
    visc_temp = visl(i,j,kp)
    dui_dxj   = (W(i,j,kp)-W(i,j,k))*dzi
    duj_dxi   = (W(i,j,kp)-W(i,j,k))*dzi
    Dz_kp     = visc_temp*(dui_dxj+duj_dxi)
    visc_temp = visl(i,j,k)
    dui_dxj   = (W(i,j,k)-W(i,j,km))*dzi
    duj_dxi   = (W(i,j,k)-W(i,j,km))*dzi
    Dz_km      = visc_temp*(dui_dxj+duj_dxi)
    dDzdz     = (Dz_kp-Dz_km)*dzi
    rhoijk    = 0.5*(rhol(i,j,k)+rhol(i,j,kp))
    diffW(i,j,k) = (1./rhoijk)*(dDxdx+dDydy+dDzdz)

    enddo
  enddo
enddo
end subroutine momdiff



end module mod_mom
