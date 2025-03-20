module mod_debug
use mod_common
use mod_param
use mod_common_mpi
use mod_common_IBM
use mod_IBM
implicit none
private
public fillVisual, PrintGhostMirror ,PrintBulkVelocities
contains
subroutine fillVisual
implicit none
integer:: i,j,k
real::ghost_tag_center,phi_ghost_center,phi_mirror_center
real:: v_tan,v_n,vv,ww,gamma_t,eps=1e-15,sgn_ny,sgn_nz,ny,nz
do k = 1,kmax
 do j = 1,jmax
  do i = 1, imax

!    v_tan = 0.
!    v_n   = 0.
!    vv    = 0.
!    ww    = 0.
!
!    sgn_ny = ny_surf_v(i,j,k)/abs(ny_surf_v(i,j,k)+eps)
!    sgn_nz = nz_surf_v(i,j,k)/abs(nz_surf_v(i,j,k)+eps)
!    if (nabs_surf_v(i,j,k).gt.1e-12) then
!      gamma_t = atan (abs(nz_surf_v(i,j,k)/(ny_surf_v(i,j,k)+eps)))
!
!       vv    = -1.
!       ww    = 0.
!       v_tan = sgn_nz*vv*sin(gamma_t)-sgn_ny*ww*cos(gamma_t)
!       v_n   = sgn_ny*vv*cos(gamma_t)+sgn_nz*ww*sin(gamma_t)
!
!      v_tan = 1.
!      v_n = 0.
!      vv   =   sgn_nz*v_tan*sin(gamma_t)+sgn_ny*v_n*cos(gamma_t)
!      ww   =  -sgn_ny*v_tan*cos(gamma_t)+sgn_nz*v_n*sin(gamma_t)
      
!    endif

     C_Po(i,j,k,1)=  PFM_phi(i,j,k)
#ifdef IBM
     Wall_IBM(i,j,k)=cell_phi_tag(i,j,k)
     C_Po(i,j,k,2)=  nx_surf(i,j,k)
     C_Po(i,j,k,3) = ny_surf(i,j,k)
     C_Po(i,j,k,4) = nz_surf(i,j,k)
     C_Po(i,j,k,5) = nabs_surf(i,j,k)
     C_Po(i,j,k,6) = cell_phi_tag(i,j,k)
#endif
  enddo
 enddo
enddo

end subroutine fillVisual



subroutine PrintGhostMirror
implicit none
real,dimension(0:itot,0:jtot,0:ktot):: x_int_all,y_int_all,z_int_all,x_mir_all,y_mir_all,z_mir_all,xxx_all,yyy_all,zzz_all
real,dimension(0:itot,0:jtot,0:ktot)::x_1_all,y_1_all,z_1_all,x_2_all,y_2_all,z_2_all
real,dimension(0:itot,0:jtot,0:ktot):: i_mir_all,j_mir_all,k_mir_all,i_all,j_all,k_all
real,dimension(0:itot,0:jtot,0:ktot)::i_1_all,j_1_all,k_1_all,i_2_all,j_2_all,k_2_all

integer:: i,j,k,ii,jj,kk,iii,jjj,kkk
integer,dimension(0:itot,0:jtot,0:ktot):: ghost_tag_all
real:: eps= 1.e-12
real:: nx,ny,nz,n_abs,x_int,y_int,z_int
real:: x_m,y_m,z_m,x_ghost,y_ghost,z_ghost
real:: xxx,yyy,zzz
#ifdef IBM
xxx_all(:,:,:)   = 0.
yyy_all(:,:,:)   = 0.
zzz_all(:,:,:)   = 0.
x_int_all(:,:,:) = 0.
y_int_all(:,:,:) = 0.
z_int_all(:,:,:) = 0.
x_mir_all(:,:,:) = 0.
y_mir_all(:,:,:) = 0.
z_mir_all(:,:,:) = 0.
ghost_tag_all(:,:,:) = 0
x_1_all(:,:,:)   = 0.
y_1_all(:,:,:)   = 0.
z_1_all(:,:,:)   = 0.
x_2_all(:,:,:)   = 0.
y_2_all(:,:,:)   = 0.
z_2_all(:,:,:)   = 0.

i_all(:,:,:)   = 0.
j_all(:,:,:)   = 0.
k_all(:,:,:)   = 0.
i_mir_all(:,:,:) = 0.
j_mir_all(:,:,:) = 0.
k_mir_all(:,:,:) = 0.
i_1_all(:,:,:)   = 0.
j_1_all(:,:,:)   = 0.
k_1_all(:,:,:)   = 0.
i_2_all(:,:,:)   = 0.
j_2_all(:,:,:)   = 0.
k_2_all(:,:,:)   = 0.


do i = 1,imax
 do j = 1,jmax
  do k = 1, kmax
   if  ( abs(nabs_surf(i,j,k)).gt.1e-12) then
    ii = i+coords(1)*imax
    jj = j+coords(2)*jmax
    kk = k
    xxx_all(ii,jj,kk)        = (ii)*dx-0.5*dx
    yyy_all(ii,jj,kk)        = (jj)*dy-0.5*dy
    ghost_tag_all(ii,jj,kk ) = 1
    zzz_all(ii,jj,kk)        = (kk     )*dz-0.5*dz
    x_int_all(ii,jj,kk)      = x_intersect(i,j,k)
    y_int_all(ii,jj,kk)      = y_intersect(i,j,k)
    z_int_all(ii,jj,kk)      = z_intersect(i,j,k)
    x_mir_all(ii,jj,kk)      = x_mirror(i,j,k)
    y_mir_all(ii,jj,kk)      = y_mirror(i,j,k)
    z_mir_all(ii,jj,kk)      = z_mirror(i,j,k)
    x_1_all(ii,jj,kk)        = x_IP1(i,j,k)
    y_1_all(ii,jj,kk)        = y_IP1(i,j,k)
    z_1_all(ii,jj,kk)        = z_IP1(i,j,k)
    x_2_all(ii,jj,kk)        = x_IP2(i,j,k)
    y_2_all(ii,jj,kk)        = y_IP2(i,j,k)
    z_2_all(ii,jj,kk)        = z_IP2(i,j,k)


    i_all(ii,jj,kk)        = ii
    j_all(ii,jj,kk)        = jj
    k_all(ii,jj,kk)        = kk
    i_mir_all(ii,jj,kk)      = i_mirror(i,j,k)*1.
    j_mir_all(ii,jj,kk)      = j_mirror(i,j,k)*1.
    k_mir_all(ii,jj,kk)      = k_mirror(i,j,k)*1.
    i_1_all(ii,jj,kk)        = i_IP1(i,j,k)*1.
    j_1_all(ii,jj,kk)        = j_IP1(i,j,k)*1.
    k_1_all(ii,jj,kk)        = k_IP1(i,j,k)*1.
    i_2_all(ii,jj,kk)        = i_IP2(i,j,k)*1.
    j_2_all(ii,jj,kk)        = j_IP2(i,j,k)*1.
    k_2_all(ii,jj,kk)        = k_IP2(i,j,k)*1.


   endif
  enddo
 enddo
enddo
call mpi_allreduce(MPI_IN_PLACE,xxx_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,yyy_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,zzz_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,x_int_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,y_int_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,z_int_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,x_mir_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,y_mir_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,z_mir_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,ghost_tag_all,(itot+1)*(jtot+1)*(ktot+1),mpi_integer,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,x_1_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,y_1_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,z_1_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,x_2_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,y_2_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,z_2_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)


call mpi_allreduce(MPI_IN_PLACE,i_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,j_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,k_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,i_mir_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,j_mir_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,k_mir_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,i_1_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,j_1_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,k_1_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,i_2_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,j_2_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(MPI_IN_PLACE,k_2_all,(itot+1)*(jtot+1)*(ktot+1),mpi_real8,mpi_sum,comm_cart,error)







if (myid.eq.0 ) then
  open(2019,file=datadir//'outputs.txt',position='append')
  i = int(imax/2)
  do j = 1,jtot
   do k = 1, ktot
     if (ghost_tag_all(i,j,k).eq.1) then
         write(2019,'(27E16.8)' ) xxx_all(i,j,k),yyy_all(i,j,k),zzz_all(i,j,k), &
                                x_int_all(i,j,k),y_int_all(i,j,k),z_int_all(i,j,k), &
                                x_mir_all(i,j,k),y_mir_all(i,j,k),z_mir_all(i,j,k), &
                                x_1_all(i,j,k),y_1_all(i,j,k),z_1_all(i,j,k), &
                                x_2_all(i,j,k),y_2_all(i,j,k),z_2_all(i,j,k), &
                                i_all(i,j,k),j_all(i,j,k),k_all(i,j,k), &
                                i_mir_all(i,j,k),j_mir_all(i,j,k),k_mir_all(i,j,k), &
                                i_1_all(i,j,k),j_1_all(i,j,k),k_1_all(i,j,k), &
                                i_2_all(i,j,k),j_2_all(i,j,k),k_2_all(i,j,k)
    endif
  enddo
 enddo
endif
close(2019)
#endif
end subroutine PrintGhostMirror


subroutine PrintBulkVelocities(time,ub,vb,wb)
implicit none
real,intent(in)::time,ub,vb,wb
open(2020,file=datadir//'Vbulk.txt',position='append')
if(myid.eq.0) write(2020,'(4E16.8)') time,ub,vb,wb
close (2020)
end subroutine PrintBulkVelocities



end module mod_debug
