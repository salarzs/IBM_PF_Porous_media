module mod_common
use mod_param
!
real ,dimension(0:i1,0:j1,0:k1) :: unew,vnew,wnew,pnew, &
                                   uo,uoo,vo,voo,wo,woo, &
                                   po,poo, phat
real ,target ,dimension(0:i1,0:j1,0:k1) :: dudt,dvdt,dwdt
real ,target ,dimension(0:i1,0:j1,0:k1) :: PFM_phi_vis
real(mytype) :: time,dt
real(mytype) ::  wi(itot+15), wj(jtot+15)
real, dimension(imax,jmax) :: xyrt
real, dimension(kmax) :: a,b,c
real :: forcextot,forceytot,forceztot
real :: u_bulk,v_bulk,w_bulk
!*********************** Phase Field Method*************
real ,dimension(0:i1,0:j1,0:k1) ::surf_tension_x, surf_tension_y,surf_tension_z
real, dimension(-5:i1+5,-5:j1+5,-5:k1+5) :: visl,rhol,PFM_phi,chem_pot,PFM_phio,PFM_phioo,psi
real ,dimension(0:jtot)::PFM_phi_agg
real::PFM_lambda
real::vel_max
real ,dimension(0:i1,0:j1) :: dPFM_bound_old,dPFM_bound_old_top
real:: alpha_PF
end module mod_common
module mod_common_IBM
use mod_param
real,allocatable::cell_u_tag(:,:,:),cell_v_tag(:,:,:),cell_w_tag(:,:,:),cell_phi_tag(:,:,:)
integer, allocatable:: Level_set(:,:,:),Level_set_u(:,:,:),Level_set_v(:,:,:),Level_set_w(:,:,:)
real, allocatable :: nx_surf(:,:,:),ny_surf(:,:,:),nz_surf(:,:,:),nabs_surf(:,:,:)
real, allocatable :: z_intersect(:,:,:),y_intersect(:,:,:), x_intersect(:,:,:)
integer,allocatable:: i_mirror(:,:,:),j_mirror(:,:,:),k_mirror(:,:,:)
real, allocatable  :: x_mirror(:,:,:), y_mirror(:,:,:), z_mirror(:,:,:)


real, allocatable :: deltan(:,:,:)
integer,allocatable :: i_IP1(:,:,:),j_IP1(:,:,:),k_IP1(:,:,:)
real, allocatable :: x_IP1(:,:,:), y_IP1(:,:,:), z_IP1(:,:,:)   !can be deallocated later
integer, allocatable :: i_IP2(:,:,:),j_IP2(:,:,:),k_IP2(:,:,:)
real, allocatable :: x_IP2(:,:,:), y_IP2(:,:,:), z_IP2(:,:,:) !can be deallocated later
real, allocatable :: WP1(:,:,:,:),WP2(:,:,:,:)
real, dimension(0:i1,0:j1,0:k1)::dPFM_boundIBM,dPFM_boundIBM_old
real:: y_12_mid,z_12_mid
real, dimension(0:i1,0:j1,0:k1):: Wall_IBM
real ,dimension(0:i1,0:j1,0:k1,1:6) :: C_Po
end module mod_common_IBM
!
module mod_common_mpi
use mpi
use decomp_2d
implicit none
integer :: myid,xhalo,yhalo,restartp,rankcw,xhalo2,yhalo2,xhalo3,yhalo3
integer :: comm_cart
!
logical periods(3),reorder
integer error,request,status(MPI_STATUS_SIZE)
integer right,rightfront,front,leftfront,left,leftback,back,rightback
integer, dimension(0:8) :: neighbor
integer, dimension(1:2) :: coords
real :: boundleftmyid,boundfrontmyid
!
end module mod_common_mpi
