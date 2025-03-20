module mod_interface
use mod_common
use mod_param
use mod_common_mpi
use mod_common_IBM
use mod_IBM
implicit none
private
public PFM_RHS,PFM_SurfaceTension,Wetting_radius,volume_of_droplet,ChemicalPotential,filtering_phi,write_time
contains
subroutine PFM_RHS(RHS_all,phi,i,j,k)
implicit none
integer,intent(in) ::i,j,k
real,intent(out)::RHS_all
real, dimension(-5:i1+5,-5:j1+5,-5:k1+5),intent(in) :: phi
real::advective
integer:: ip,jp,kp,im,jm,km
real::phi_p,phi_m
real::advectivex,advectivey, advectivez
real:: Pot_ijk,Pot_ip,Pot_jp,Pot_kp,Pot_im,Pot_jm,Pot_km
real:: Cahn_Hilliard_RHS,d2potdx2,d2potdy2,d2potdz2

ip=i+1
jp=j+1
kp=k+1
im=i-1
jm=j-1
km=k-1

advective=0.0
Cahn_Hilliard_RHS = 0.
RHS_all=0.0
advectivex=0.0
advectivey=0.0
advectivez=0.0
phi_p= 0.0
phi_m= 0.0

phi_p=phi(ip,j,k)
phi_m=phi(im,j,k)
advectivex=-0.5*(unew(i,j,k)+unew(i-1,j,k))*(phi_p-phi_m)*(0.5*dxi)

phi_p= 0.0
phi_m= 0.0
phi_p=phi(i,jp,k)
phi_m=phi(i,jm,k)
advectivey=-0.5*(vnew(i,j,k)+vnew(i,j-1,k))*(phi_p-phi_m)*(0.5*dyi)

phi_p= 0.0
phi_m= 0.0
phi_p=phi(i,j,kp)
phi_m=phi(i,j,km)
advectivez=-0.5*(wnew(i,j,k)+wnew(i,j,k-1))*(phi_p-phi_m)*(0.5*dzi)

advective=advectivex+advectivey+advectivez

Pot_ijk  = chem_pot(i,j,k)
Pot_ip   = chem_pot(ip,j,k)
Pot_im   = chem_pot(im,j,k)
Pot_jp   = chem_pot(i,jp,k)
Pot_jm   = chem_pot(i,jm,k)
Pot_kp   = chem_pot(i,j,kp)
Pot_km   = chem_pot(i,j,km)
d2potdx2 = (Pot_ip-2.*Pot_ijk+Pot_im)*dxi*dxi
d2potdy2 = (Pot_jp-2.*Pot_ijk+Pot_jm)*dyi*dyi
d2potdz2 = (Pot_kp-2.*Pot_ijk+Pot_km)*dzi*dzi
if (Phi01)      Cahn_Hilliard_RHS =2* mobility*(d2potdx2+d2potdy2+d2potdz2)    ! For 0<phi<1
if (.not.phi01) Cahn_Hilliard_RHS = mobility*(d2potdx2+d2potdy2+d2potdz2)
RHS_all = Cahn_Hilliard_RHS  +advective
#ifdef IBM
    if (nabs_surf(i,j,k).gt.1e-12) RHS_all = 0.
#endif
return
end subroutine PFM_RHS


subroutine PFM_grad_phi(DphiDx,DphiDy,DphiDz,Phi_grad_abs,i,j,k,mode)
implicit none
integer,intent(in) ::i,j,k,mode
real,intent(out)::DphiDx
real,intent(out)::DphiDy
real,intent(out)::DphiDz
real,intent(out)::Phi_grad_abs
real:: phi_ip,phi_im,phi_jp,phi_jm ,phi_kp,phi_km,phi_kpp,phi_kmm,phi_ijk

DphiDx=0.0
DphiDy=0.0
DphiDz=0.0
Phi_grad_abs=0.0
phi_ip=0.0
phi_im=0.0
phi_jp=0.0
phi_jm=0.0
phi_kp=0.0
phi_km =0.0
phi_kpp=0.0
phi_kmm=0.0

! at the center
if (mode.eq.0) then
    DphiDx=(PFM_phi(i+1,j,k)-PFM_phi(i-1,j,k))*(0.5*dxi)
    DphiDy=(PFM_phi(i,j+1,k)-PFM_phi(i,j-1,k))*(0.5*dyi)
    DphiDz=(PFM_phi(i,j,k+1)-PFM_phi(i,j,k-1))*(0.5*dzi)
 
    Phi_grad_abs=sqrt(DphiDx*DphiDx+DphiDy*DphiDy+DphiDz*DphiDz)
    if (Phi_grad_abs.ne.0) then
    DphiDx=DphiDx/Phi_grad_abs
    DphiDy=DphiDy/Phi_grad_abs
    DphiDz=DphiDz/Phi_grad_abs
    endif
    
    if (Phi_grad_abs.eq.0) then
    DphiDx= 0.0
    DphiDy= 0.0
    DphiDz= 0.0
    endif
endif
! at u face

if (mode.eq.1) then
    phi_ijk  = 0.5*(PFM_phi(i+1,j,k)+PFM_phi(i,j,k))
    phi_ip   = 0.5*(PFM_phi(i+2,j,k)+PFM_phi(i+1,j,k))
    phi_im   = 0.5*(PFM_phi(i,j,k)+PFM_phi(i-1,j,k))
    phi_jp   = 0.5*(PFM_phi(i+1,j+1,k)+PFM_phi(i,j+1,k))
    phi_jm   = 0.5*(PFM_phi(i+1,j-1,k)+PFM_phi(i,j-1,k))
    phi_kp   = 0.5*(PFM_phi(i+1,j,k+1)+PFM_phi(i,j,k+1))
    phi_km   = 0.5*(PFM_phi(i+1,j,k-1)+PFM_phi(i,j,k-1))


    DphiDx=(phi_ip-phi_im)*(0.5*dxi)
    DphiDy=(phi_jp-phi_jm)*(0.5*dyi)
    DphiDz=(phi_kp-phi_km)*(0.5*dzi) 

    Phi_grad_abs=sqrt(DphiDx*DphiDx+DphiDy*DphiDy+DphiDz*DphiDz)
    if (Phi_grad_abs.ne.0) then
    DphiDx=DphiDx/Phi_grad_abs
    DphiDy=DphiDy/Phi_grad_abs
    DphiDz=DphiDz/Phi_grad_abs
    endif
    
    if (Phi_grad_abs.eq.0) then
    DphiDx= 0.0
    DphiDy= 0.0
    DphiDz= 0.0
    endif

endif

! at v face

if (mode.eq.2) then
    phi_ip = 0.5*(PFM_phi(i+1,j+1,k)+PFM_phi(i+1,j,k))
    phi_im = 0.5*(PFM_phi(i-1,j+1,k)+PFM_phi(i-1,j,k))
    phi_jp = 0.5*(PFM_phi(i,j+2,k)+PFM_phi(i,j+1,k))
    phi_jm = 0.5*(PFM_phi(i,j,k)+PFM_phi(i,j-1,k))
    phi_kp = 0.5*(PFM_phi(i,j+1,k+1)+PFM_phi(i,j,k+1))
    phi_km = 0.5*(PFM_phi(i,j+1,k-1)+PFM_phi(i,j,k-1))
    
    DphiDx=(phi_ip-phi_im)*(0.5*dxi)
    DphiDy=(phi_jp-phi_jm)*(0.5*dyi)
    DphiDz=(phi_kp-phi_km)*(0.5*dzi)

    Phi_grad_abs=sqrt(DphiDx*DphiDx+DphiDy*DphiDy+DphiDz*DphiDz)
    if (Phi_grad_abs.ne.0) then
    DphiDx=DphiDx/Phi_grad_abs
    DphiDy=DphiDy/Phi_grad_abs
    DphiDz=DphiDz/Phi_grad_abs
    endif
    
    if (Phi_grad_abs.eq.0) then
    DphiDx= 0.0
    DphiDy= 0.0
    DphiDz= 0.0
    endif
    
endif
 ! at w face

if (mode.eq.3) then
    phi_ip = 0.5*(PFM_phi(i+1,j,k+1)+PFM_phi(i+1,j,k))
    phi_im = 0.5*(PFM_phi(i-1,j,k+1)+PFM_phi(i-1,j,k))
    phi_jp = 0.5*(PFM_phi(i,j+1,k+1)+PFM_phi(i,j+1,k))
    phi_jm = 0.5*(PFM_phi(i,j-1,k+1)+PFM_phi(i,j-1,k))
    phi_kp = 0.5*(PFM_phi(i,j,k+2)+PFM_phi(i,j,k+1))
    phi_km = 0.5*(PFM_phi(i,j,k)+PFM_phi(i,j,k-1))
    
    DphiDx=(phi_ip-phi_im)*(0.5*dxi)
    DphiDy=(phi_jp-phi_jm)*(0.5*dyi)
    DphiDz=(phi_kp-phi_km)*(0.5*dzi)

    Phi_grad_abs=sqrt(DphiDx*DphiDx+DphiDy*DphiDy+DphiDz*DphiDz)
    if (Phi_grad_abs.ne.0) then
    DphiDx=DphiDx/Phi_grad_abs
    DphiDy=DphiDy/Phi_grad_abs
    DphiDz=DphiDz/Phi_grad_abs
    endif
    
    if (Phi_grad_abs.eq.0) then
    DphiDx= 0.0
    DphiDy= 0.0
    DphiDz= 0.0
    endif
    
endif
return
end subroutine PFM_grad_phi


subroutine PFM_SurfaceTension (i,j,k,Fx_SurfTen,Fy_SurfTen,Fz_SurfTen)
implicit none
integer,intent(in) ::i,j,k
real,intent(out)::Fx_SurfTen,Fy_SurfTen,Fz_SurfTen
integer:: ip,jp,kp,im,jm,km
real::Pot_ijk,Pot_ip,Pot_jp,Pot_kp
real:: rhoi,rhoj,rhok,phii,phij,phik
real::DphiDx,DphiDy,DphiDz,Phi_grad_abs

ip=i+1
jp=j+1
kp=k+1
im=i-1
jm=j-1
km=k-1

Fx_SurfTen=0.0
Fy_SurfTen=0.0
Fz_SurfTen=0.0


phii = 0.5*(PFM_phi(i,j,k)+PFM_phi(ip,j,k))
phij = 0.5*(PFM_phi(i,j,k)+PFM_phi(i,jp,k))
phik = 0.5*(PFM_phi(i,j,k)+PFM_phi(i,j,kp))


rhoi= 0.5*(phii+1.)*rho_2-0.5*(phii-1.)*rho_1
rhoj= 0.5*(phij+1.)*rho_2-0.5*(phij-1.)*rho_1
rhok= 0.5*(phik+1.)*rho_2-0.5*(phik-1.)*rho_1

Pot_ijk = chem_pot(i,j,k)
Pot_ip  = chem_pot(ip,j,k)
Pot_jp  = chem_pot(i,jp,k)
Pot_kp  = chem_pot(i,j,kp)

call PFM_grad_phi(DphiDx,DphiDy,DphiDz,Phi_grad_abs,i,j,k,1)
Fx_SurfTen = 0.5*(Pot_ijk+Pot_ip)*DphiDx*Phi_grad_abs/rhoi
call PFM_grad_phi(DphiDx,DphiDy,DphiDz,Phi_grad_abs,i,j,k,2)
Fy_SurfTen = 0.5*(Pot_ijk+Pot_jp)*DphiDy*Phi_grad_abs/rhoj
call PFM_grad_phi(DphiDx,DphiDy,DphiDz,Phi_grad_abs,i,j,k,3)
Fz_SurfTen = 0.5*(Pot_ijk+Pot_kp)*DphiDz*Phi_grad_abs/rhok
 
return

end subroutine PFM_SurfaceTension

subroutine ChemicalPotential
implicit none
integer::i,j,k
integer:: ip,im,jp,jm,kp,km
real::phi_ijk,phi_ip,phi_jp,phi_kp,phi_im,phi_jm,phi_km
real::Psi_prime,d2phidx2,d2phidy2,d2phidz2,lap_phi

do k=1,kmax
 do j=1,jmax
   do i=1,imax
#ifdef IBM
    if (abs(nabs_surf(i,j,k)).gt.1e-12) cycle
#endif
     ip=i+1
     im=i-1
     jp=j+1
     jm=j-1
     kp=k+1
     km=k-1
     phi_ijk = PFM_phi(i,j,k)
     phi_ip  = PFM_phi(ip,j,k)
     phi_jp  = PFM_phi(i,jp,k)
     phi_kp  = PFM_phi(i,j,kp)
     phi_im  = PFM_phi(im,j,k)
     phi_jm  = PFM_phi(i,jm,k)
     phi_km  = PFM_phi(i,j,km)
     if (.not.Phi01)   Psi_prime = phi_ijk*(phi_ijk*phi_ijk-1.)     ! For -1<phi<1
     if (Phi01)        Psi_prime =phi_ijk*(phi_ijk*phi_ijk-1.)/8. ! For 0<phi<1
     d2phidx2 = (phi_ip-2.*phi_ijk+phi_im)*dxi*dxi
     d2phidy2 = (phi_jp-2.*phi_ijk+phi_jm)*dyi*dyi
     d2phidz2 = (phi_kp-2.*phi_ijk+phi_km)*dzi*dzi
     lap_phi= d2phidx2 + d2phidy2 + d2phidz2
     if (.not.Phi01) Chem_Pot(i,j,k) =(3./(2.*sqrt(2.)))*((PFM_Sigma_ref/PFM_l)*(Psi_prime)-(PFM_Sigma_ref*PFM_l*lap_phi)) ! For -1<phi<1
     if (Phi01) Chem_Pot(i,j,k) =(3./(2.*sqrt(2.)))*((PFM_Sigma_ref/PFM_l)*(Psi_prime)-(PFM_Sigma_ref*PFM_l*lap_phi*0.5)) ! For 0<phi<1
   enddo
  enddo
enddo
return

end subroutine ChemicalPotential

subroutine Wetting_radius(time)

implicit none
integer::k,j,i
real:: starting_point_j,end_point_j,starting_point_k,end_point_k,phi_wall
real,intent(in)::time
real:: time_stari_visc,time_star_inert
real::Wetting_rad, phi_m,delta_y,delta_z 
logical:: confirmation
real,dimension(1:dims(2)):: j_min,j_max,k_min,k_max
integer::min_processor,max_processor
real::j_min_all,j_max_all,dn_1,dn_2,alpha
time_star_inert= sqrt((rho_2*droplet_radius*droplet_radius*droplet_radius)/PFM_sigma_ref)
time_star_inert= time/time_star_inert

time_stari_visc = vis_2*droplet_radius/PFM_sigma_ref
time_stari_visc = time/time_stari_visc

j_min(:) =  10000.
j_max(:) = -10000.
k_min(:) =  10000.
k_max(:) = -10000.
Wetting_rad = 0.
j_min_all =  10000.
j_max_all = -10000.
min_processor = 0
max_processor = 0
#ifdef NIBM
   i= imax/2
   do j=1,jmax
         phi_wall = 0.5*(PFM_phi(i,j,0)+PFM_phi(i,j,1))
         if ((phi_wall).gt.0.0) then
         j_min(coords(2)+1)= (j+coords(2)*jmax)*1.0
         k_min(coords(2)+1)= k*1.0
         exit
         endif
   enddo
   
   do j=jmax,1,-1
         phi_wall = 0.5*(PFM_phi(i,j,0)+PFM_phi(i,j,1))
         if ((phi_wall).gt.0.0) then
         j_max(coords(2)+1)=  (j+coords(2)*jmax)*1.0
         k_max(coords(2)+1)= k *1.0 
         exit
         endif
   enddo

#endif
#ifdef IBM
  confirmation = .false.
  delta_y = 0.
  delta_z = 0.
  j_min(:) =  10000.
  j_max(:) = -10000.
  k_min(:) =  10000.
  k_max(:) = -10000.

  Wetting_rad = 0.
  j_min_all = 0.
  j_max_all = 0.
  min_processor = 0
  max_processor = 0
  i=int(imax/2)
  do j=1,jmax
    do k=1,kmax
      if  ( abs(nabs_surf(i,j,k)).gt.1e-12) then
          call interpolation_mirror(PFM_phi,i,j,k,phi_m)
          if ((phi_m).gt.0.0) then
           j_min(coords(2)+1)=  (j+coords(2)*jmax)*1.0
           k_min(coords(2)+1)= k*1.0 
           confirmation = .true.
           exit 
          endif
       endif
    enddo
    if (confirmation) exit
  enddo
  
  confirmation = .false.
  do j=jmax,1,-1
   do k=kmax,1,-1
      if  ( abs(nabs_surf(i,j,k)).gt.1e-12) then
          call interpolation_mirror(PFM_phi,i,j,k,phi_m)
          if ((phi_m).gt.0.0) then
           j_max(coords(2)+1)= (j+coords(2)*jmax)*1.0
           k_max(coords(2)+1)= k*1.0
           confirmation = .true.
           exit
          endif
      endif
   enddo
   if (confirmation) exit
  enddo

  call mpi_allreduce(MPI_IN_PLACE,j_min(1),dims(2),mpi_real8,mpi_min,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,j_max(1),dims(2),mpi_real8,mpi_max,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,k_min(1),dims(2),mpi_real8,mpi_min,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,k_max(1),dims(2),mpi_real8,mpi_max,comm_cart,error)



  
j_min_all =  10000.
j_max_all = -10000.
  do i = 1, dims(2)
    if (j_min(i).lt.j_min_all) then
        j_min_all = j_min(i)
        min_processor = i
    endif
    if (j_max(i).gt.j_max_all) then
        j_max_all = j_max(i)
        max_processor = i
    endif
enddo

  delta_y = (j_max_all-j_min_all)*dy
  delta_z = (k_max(max_processor)-k_min(min_processor))*dz
  Wetting_rad = 0.5*sqrt(delta_y*delta_y+delta_z*delta_z)/droplet_radius
  open(1329,file=datadir//'spreading_radius.txt',position='append')
  if (myid.eq.0)  write(1329,'(3E16.8)' ) time_stari_visc,time_star_inert,Wetting_rad
  close(1329)






  confirmation = .false.
  delta_y = 0.
  delta_z = 0.
  j_min(:) =  10000.
  j_max(:) = -10000.
  k_min(:) =  10000.
  k_max(:) = -10000.

  Wetting_rad = 0.
  j_min_all = 0.
  j_max_all = 0.
  min_processor = 0
  max_processor = 0
  i=int(imax/2)
  do j=1,jmax
    do k=1,kmax
      if  ( abs(nabs_surf(i,j,k)).gt.1e-12) then
          call interpolation_mirror(PFM_phi,i,j,k,phi_m)
          dn_1 = deltan(i,j,k)
          dn_2 = march_step
          alpha= dn_2/dn_1
          phi_wall = (alpha*PFM_phi(i,j,k)+phi_m)/(alpha+1)
          if ((phi_wall).gt.0.0) then
           j_min(coords(2)+1)=  (j+coords(2)*jmax)*1.0
           k_min(coords(2)+1)= k*1.0
           confirmation = .true.
           exit
          endif
       endif
    enddo
    if (confirmation) exit
  enddo

  confirmation = .false.
  do j=jmax,1,-1
   do k=kmax,1,-1
      if  ( abs(nabs_surf(i,j,k)).gt.1e-12) then
          call interpolation_mirror(PFM_phi,i,j,k,phi_m)
          dn_1 = deltan(i,j,k)
          dn_2 = march_step
          alpha= dn_2/dn_1
          phi_wall = (alpha*PFM_phi(i,j,k)+phi_m)/(alpha+1)
          if ((phi_wall).gt.0.0) then
           j_max(coords(2)+1)= (j+coords(2)*jmax)*1.0
           k_max(coords(2)+1)= k*1.0
           confirmation = .true.
           exit
          endif
      endif
   enddo
   if (confirmation) exit
  enddo

  call mpi_allreduce(MPI_IN_PLACE,j_min(1),dims(2),mpi_real8,mpi_min,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,j_max(1),dims(2),mpi_real8,mpi_max,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,k_min(1),dims(2),mpi_real8,mpi_min,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,k_max(1),dims(2),mpi_real8,mpi_max,comm_cart,error)




j_min_all =  10000.
j_max_all = -10000.
  do i = 1, dims(2)
    if (j_min(i).lt.j_min_all) then
        j_min_all = j_min(i)
        min_processor = i
    endif
    if (j_max(i).gt.j_max_all) then
        j_max_all = j_max(i)
        max_processor = i
    endif
enddo

  delta_y = (j_max_all-j_min_all)*dy
  delta_z = (k_max(max_processor)-k_min(min_processor))*dz
  Wetting_rad = 0.5*sqrt(delta_y*delta_y+delta_z*delta_z)/droplet_radius
  open(1328,file=datadir//'spreading_radius2.txt',position='append')
  if (myid.eq.0)  write(1328,'(3E16.8)' ) time_stari_visc,time_star_inert,Wetting_rad
  close(1328)
#endif
return
end subroutine Wetting_radius

subroutine write_time(istep,time)
implicit none
real,intent(in)::time
integer,intent(in):: istep
real:: time_star
time_star = time*u_ref/droplet_radius
open(1327,file=datadir//'time_step',position='append')
if (myid.eq.0)  write(1327,'(3E16.8)' ) time, time_star, 1.*istep
close(1327)
return
end subroutine write_time


subroutine volume_of_droplet(time)
implicit none
integer:: i,j,k
real:: v_droplet 
real,intent(in)::time

v_droplet = 0.
do k=1,kmax
 do j=1,jmax
   do i=1,imax
#ifdef NIBM
        v_droplet= v_droplet+ dx*dy*dz*0.5*(PFM_phi(i,j,k)+1)
#endif
#ifdef IBM
        v_droplet= v_droplet+ dx*dy*dz*0.5*(PFM_phi(i,j,k)+1)*cell_phi_tag(i,j,k)
#endif
   enddo
 enddo
enddo
call mpi_allreduce(MPI_IN_PLACE,v_droplet,1,mpi_real8,mpi_sum,comm_cart,error)
open(1341,file=datadir//'volume_of_droplet.txt',position='append')
if (myid.eq.0)  write(1341,'(2E16.8)' ) time,v_droplet
close(1341)
return
end subroutine volume_of_droplet


subroutine filtering_phi
implicit none
integer:: i,j,k
  do k=-5,k1+5
   do j=-5,j1+5
    do i=-5,i1+5
      PFM_phi(i,j,k)=MIN(MAX(PFM_phi(i,j,k),-1.0),1.0)
      rhol(i,j,k)= 0.5*(PFM_phi(i,j,k)+1)*rho_2-0.5*(PFM_phi(i,j,k)-1)*rho_1
      visl(i,j,k)= 0.5*(PFM_phi(i,j,k)+1)*vis_2-0.5*(PFM_phi(i,j,k)-1)*vis_1
    enddo
   enddo
  enddo
end subroutine filtering_phi
end module mod_interface 
