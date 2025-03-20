module mod_init
use mod_common
use decomp_2d
use mod_param
use mod_common_mpi
use mod_common_IBM
use mod_IBM
use mod_bound
implicit none
private
public init
contains
subroutine init
implicit none
real :: distance,xxx,yyy,zzz,DropletsCentx, DropletsCenty,DropletsCentz
integer::i,j,k
real::phi_temp,y_p,z_p,y_prime,z_prime,y_tilde,z_tilde
real,dimension(0:i1,0:j1,0:k1):: v_trans,w_trans,vcenter,wcenter
real::d
unew(:,:,:)    = 0.
vnew(:,:,:)    = 0.
wnew(:,:,:)    = 0.
pnew(:,:,:)    = 0.
PFM_phi(:,:,:) = 0.
select case(PhaseField)
case('Droplet')
PFM_phi(:,:,:)=-1
#ifdef TwoD
#ifdef IBM
DropletsCenty= ly/2. !ly/2 !(droplet_radius+0.5*(Droplets_separation-2*droplet_radius))
DropletsCentz=droplet_radius + (h_w)*dz+0.2*lz
do k=0,k1
  do j=0,j1
   do i=0,i1
      xxx=(i+coords(1)*imax-1)*dx+0.5*dx
      yyy=(j+coords(2)*jmax-1)*dy+0.5*dy
      zzz=(k               -1)*dz+0.5*dz
      PFM_phi(i,j,k)=-PFM_phi_range
      distance= 0.
      distance= sqrt(( yyy-DropletsCenty )**2+(zzz-DropletsCentz )**2)
      PFM_phi(i,j,k) = -PFM_phi_range*tanh((distance-droplet_radius)/(PFM_l*sqrt(2.)))
      if (IBM) then
      if ((cell_phi_tag(i,j,k).eq.0).and.(abs(nabs_surf(i,j,k)).lt.1e-12)) then
          PFM_phi(i,j,k)=-1
      endif

      endif
   enddo
  enddo
enddo



#endif
#ifdef NIBM
DropletsCenty= ly/2. 
DropletsCentz=  0.95*droplet_radius
do k=0,k1
  do j=0,j1
   do i=0,i1
      xxx=(i+coords(1)*imax-1)*dx+0.5*dx
      yyy=(j+coords(2)*jmax-1)*dy+0.5*dy
      zzz=(k               -1)*dz+0.5*dz
      PFM_phi(i,j,k)=-PFM_phi_range
      distance= 0.
      distance= sqrt(( yyy-DropletsCenty )**2+(zzz-DropletsCentz )**2)
      PFM_phi(i,j,k) = -PFM_phi_range*tanh((distance-droplet_radius)/(PFM_l*sqrt(2.)))
   enddo
  enddo
enddo
#endif

#else
#ifdef IBM
DropletsCentx= lx/2.
DropletsCenty= ly/2. !ly/2 !(droplet_radius+0.5*(Droplets_separation-2*droplet_radius))
DropletsCentz=  droplet_radius+(h_w+1.*h_b)*dz

do k=0,k1
  do j=0,j1
   do i=0,i1
      xxx=(i+coords(1)*imax-1)*dx+0.5*dx
      yyy=(j+coords(2)*jmax-1)*dy+0.5*dy
      zzz=(k               -1)*dz+0.5*dz
      PFM_phi(i,j,k)=-PFM_phi_range
      distance= 0.
      distance= sqrt((xxx-DropletsCentx)**2.+( yyy-DropletsCenty)**2.+(zzz-DropletsCentz )**2.)
      PFM_phi(i,j,k) = -PFM_phi_range*tanh((distance-droplet_radius)/(PFM_l*sqrt(2.)))
       if ((cell_phi_tag(i,j,k).eq.0).and.(abs(nabs_surf(i,j,k)).lt.1e-12)) then
           PFM_phi(i,j,k)=-1
      endif

         if(yyy.gt.0.5*ly) then
            thetta_array(i,j) = 45.*pi/180
         else
            thetta_array(i,j) = 135.*pi/180
         endif 

   enddo
  enddo
enddo


#endif
#ifdef NIBM
DropletsCentx= lx/2.
DropletsCenty= ly/2. !ly/2 !(droplet_radius+0.5*(Droplets_separation-2*droplet_radius))
DropletsCentz= 0.9*droplet_radius 
do k=0,k1
  do j=0,j1
   do i=0,i1
      xxx=(i+coords(1)*imax-1)*dx+0.5*dx
      yyy=(j+coords(2)*jmax-1)*dy+0.5*dy
      zzz=(k               -1)*dz+0.5*dz
      PFM_phi(i,j,k)=-PFM_phi_range
       distance= 0.
      distance= sqrt((xxx-DropletsCentx)**2.+( yyy-DropletsCenty)**2.+(zzz-DropletsCentz )**2.)
      PFM_phi(i,j,k) = -PFM_phi_range*tanh((distance-droplet_radius)/(PFM_l*sqrt(2.)))
   enddo
  enddo
enddo

#endif

#endif

 
call updthalosBig(PFM_phi,1)
call updthalosBig(PFM_phi,2)
end select
end subroutine init
!
end module mod_init
