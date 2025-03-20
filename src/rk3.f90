module mod_rk
use mod_mom
use mod_common
use mod_common_mpi
use mod_bound
use mod_common_IBM
use mod_IBM
implicit none
private
public updatePfm,ab2,ab2_Imp,updatePfmDyn
contains
!
subroutine updatePfm(dPFM)
implicit none
integer i,j,k,ii
real :: factor1,factor2
real ::dPFMnew(0:i1,0:j1,0:k1)
real, intent(inout) :: dPFM(0:i1,0:j1,0:k1)
real:: Phi_wall_n(0:i1,0:j1),phi_wall_n_top(0:i1,0:j1)
real::phi_wall_nIBM(0:i1,0:j1,0:k1)
real:: dPFM_bound_new,Phi_wall_n1,Phi_wall_n1_top
real, dimension(0:i1,0:j1) ::dPFM_bound_n,dPFM_bound_n_top
real:: nx,ny,nz,n_abs,phi_m
integer:: i_p, j_p, k_p
real:: eps = 1e-12
integer:: i_mir,j_mir,k_mir
real:: xxx,yyy,zzz
real::dn_1,dn_2,alpha
factor1 = dt*( 1.5)
factor2 = dt*(-0.5)

#ifdef NIBM
  Phi_wall_n=0.
  phi_wall_n_top=0.
  !at time n:
  do j=1,jmax
    do i=1,imax
      phi_wall_n(i,j)=0.5*(PFM_phi(i,j,1)+PFM_phi(i,j,0))
      phi_wall_n_top(i,j)=0.5*(PFM_phi(i,j,kmax)+PFM_phi(i,j,k1))
    enddo
  enddo
  call boundPFM(PFM_phi,dPFM_bound_n,dPFM_bound_n_top,unew,vnew,wnew)
  !update phi to time n+1
  call PFM(dPFMnew,PFM_phi)
  do k=1,kmax
    do j=1,jmax
      do i=1,imax
        PFM_phi(i,j,k) =PFM_phi(i,j,k)+factor2*dPFM(i,j,k) + factor1*dPFMnew(i,j,k)
        dPFM(i,j,k)=dPFMnew(i,j,k)
      enddo
    enddo
  enddo
  !update the boundary to n+1
  do j=1,jmax
    do i=1,imax
      Phi_wall_n1=phi_wall_n(i,j)+factor1*dPFM_bound_n(i,j)+factor2*dPFM_bound_old(i,j)
      Phi_wall_n1_top=phi_wall_n_top(i,j)+factor1*dPFM_bound_n_top(i,j)+factor2*dPFM_bound_old_top(i,j)
      if (abs(PFM_mu_f).gt.1e-12) then
       PFM_phi(i,j,0) = 2*Phi_wall_n1-PFM_phi(i,j,1)
       dPFM_bound_old(i,j)=dPFM_bound_n(i,j)
       select case(PhaseField)
        case ('Couette')
          PFM_phi(i,j,k1) = 2*Phi_wall_n1_top-PFM_phi(i,j,kmax)
          dPFM_bound_old_top(i,j)=dPFM_bound_n_top(i,j)
       end select
      endif
    enddo
  enddo
#endif

#ifdef IBM
   do k=1,kmax
    do j=1,jmax
     do i=1,imax
         phi_wall_nIBM(i,j,k) = 0.
         nx =  nx_surf(i,j,k)
         ny =  ny_surf(i,j,k)
         nz =  nz_surf(i,j,k)
         n_abs = nabs_surf(i,j,k)
         dn_1 = deltan(i,j,k)
         dn_2 = march_step
         alpha= dn_2/(dn_1+1e-12)

         if  ( abs(nabs_surf(i,j,k)).gt.1e-12) then
               call interpolation_mirror(PFM_phi,i,j,k,phi_m)
               xxx=(i+coords(1)*imax)*dx-0.5*dx
               yyy=(j+coords(2)*jmax)*dy-0.5*dy
               zzz=(k               )*dz-0.5*dz
               phi_wall_nIBM(i,j,k) = (alpha*PFM_phi(i,j,k)+phi_m)/(alpha+1)
         endif
     enddo
    enddo
   enddo
   call boundPFM(PFM_phi,dPFM_bound_n,dPFM_bound_n_top,unew,vnew,wnew)   
   !update phi to time n+1
   call PFM(dPFMnew,PFM_phi)
   do k=1,kmax
     do j=1,jmax
       do i=1,imax
         PFM_phi(i,j,k) =PFM_phi(i,j,k)+factor2*dPFM(i,j,k) + factor1*dPFMnew(i,j,k)
         dPFM(i,j,k)=dPFMnew(i,j,k)
       enddo
      enddo
  enddo
   !update the boundary to n+1
   do k=1,kmax
    do j=1,jmax
     do i=1,imax
         ! Inside the wall
         if  ( abs(nabs_surf(i,j,k)).gt.1e-12) then
           nx =  nx_surf(i,j,k)
           ny =  ny_surf(i,j,k)
           nz =  nz_surf(i,j,k)
           n_abs = sqrt(nx*nx+ny*ny+nz*nz)
           Phi_wall_n1=phi_wall_nIBM(i,j,k) +factor1*dPFM_boundIBM(i,j,k)+factor2*dPFM_boundIBM_old(i,j,k)
   
           if (abs(PFM_mu_f).gt.1e-12) then
              xxx=(i+coords(1)*imax)*dx-0.5*dx
              yyy=(j+coords(2)*jmax)*dy-0.5*dy
              zzz=(k               )*dz-0.5*dz
              call interpolation_mirror(PFM_phi,i,j,k,phi_m)
              PFM_phi(i,j,k) = ((alpha+1.)*Phi_wall_n1-phi_m)/(alpha)
              dPFM_boundIBM_old(i,j,k)=dPFM_boundIBM(i,j,k)
           endif
         endif
     enddo
    enddo
   enddo

#endif

end subroutine updatePfm
!
subroutine ab2(durk1,dvrk1,dwrk1,dPFM_couple_x,dPFM_couple_y,dPFM_couple_z)
implicit none
integer i,j,k,ii
real, intent(inout) :: durk1(0:,0:,0:)
real, intent(inout) :: dvrk1(0:,0:,0:)
real, intent(inout) :: dwrk1(0:,0:,0:)
real :: dnewu(0:i1,0:j1,0:k1) ! dummy array
real :: dnewv(0:i1,0:j1,0:k1) ! dummy array
real :: dneww(0:i1,0:j1,0:k1) ! dummy array
real :: durkDiff(0:i1,0:j1,0:k1) ! dummy array
real :: dvrkDiff(0:i1,0:j1,0:k1) ! dummy array
real :: dwrkDiff(0:i1,0:j1,0:k1) ! dummy array

real ::dPFM_couple_xnew(0:i1,0:j1,0:k1)
real ::dPFM_couple_ynew(0:i1,0:j1,0:k1)
real ::dPFM_couple_znew(0:i1,0:j1,0:k1)
real, intent(inout) :: dPFM_couple_x(0:,0:,0:)
real, intent(inout) :: dPFM_couple_y(0:,0:,0:)
real, intent(inout) :: dPFM_couple_z(0:,0:,0:)
real:: duDiffVof(0:i1,0:j1,0:k1)
real:: dvDiffVof(0:i1,0:j1,0:k1)
real:: dwDiffVof(0:i1,0:j1,0:k1)
real :: factor1,factor2
!----------- Two-dimensional flow--------
factor1 = dt*( 1.5)
factor2 = dt*(-0.5)

do k=0,k1
  do j=0,j1
    do i=0,i1
      duDiffVof(i,j,k) = 0.
      dvDiffVof(i,j,k) = 0.
      dwDiffVof(i,j,k) = 0.

    enddo
  enddo
enddo

call momxad(dnewu,durkDiff,unew,vnew,wnew)
call momyad(dnewv,dvrkDiff,unew,vnew,wnew)
call momzad(dneww,dwrkDiff,unew,vnew,wnew)
call momdiffvof(duDiffVof,dvDiffVof,dwDiffVof,unew,vnew,wnew,PFM_phi)
do k=1,kmax
  do j=1,jmax
    do i=1,imax
     durkDiff(i,j,k)=durkDiff(i,j,k)+duDiffVof(i,j,k)
     dvrkDiff(i,j,k)=dvrkDiff(i,j,k)+dvDiffVof(i,j,k)
     dwrkDiff(i,j,k)=dwrkDiff(i,j,k)+dwDiffVof(i,j,k)


      dnewu(i,j,k)= dnewu(i,j,k)+durkDiff(i,j,k)
      dnewv(i,j,k)= dnewv(i,j,k)+dvrkDiff(i,j,k)
      dneww(i,j,k)= dneww(i,j,k)+dwrkDiff(i,j,k) 

      dudt(i,j,k) = unew(i,j,k) + factor1*dnewu(i,j,k) + factor2*durk1(i,j,k) 
      dvdt(i,j,k) = vnew(i,j,k) + factor1*dnewv(i,j,k) + factor2*dvrk1(i,j,k) 
      dwdt(i,j,k) = wnew(i,j,k) + factor1*dneww(i,j,k) + factor2*dwrk1(i,j,k) 

      durk1(i,j,k) = dnewu(i,j,k)
      dvrk1(i,j,k) = dnewv(i,j,k)
      dwrk1(i,j,k) = dneww(i,j,k)
    enddo
  enddo
enddo

  call PFM_couple_x(dPFM_couple_xnew)
  call PFM_couple_y(dPFM_couple_ynew)
  call PFM_couple_z(dPFM_couple_znew)
  do k=1,kmax
    do j=1,jmax
      do i=1,imax
        dudt(i,j,k) = dudt(i,j,k) + factor1*dPFM_couple_xnew(i,j,k)+ factor2*dPFM_couple_x(i,j,k)
        dvdt(i,j,k) = dvdt(i,j,k) + factor1*dPFM_couple_ynew(i,j,k)+ factor2*dPFM_couple_y(i,j,k)
        dwdt(i,j,k) = dwdt(i,j,k) + factor1*dPFM_couple_znew(i,j,k)+ factor2*dPFM_couple_z(i,j,k)
       
        dPFM_couple_x(i,j,k)=dPFM_couple_xnew(i,j,k)
        dPFM_couple_y(i,j,k)=dPFM_couple_ynew(i,j,k)
        dPFM_couple_z(i,j,k)=dPFM_couple_znew(i,j,k)

        surf_tension_x(i,j,k)= dPFM_couple_x(i,j,k)
        surf_tension_y(i,j,k)= dPFM_couple_y(i,j,k)
        surf_tension_z(i,j,k)= dPFM_couple_z(i,j,k)
        
       
      enddo
    enddo
  enddo
#ifdef TwoD
dudt(:,:,:) = 0.
unew(:,:,:) = 0.
#endif
return
end subroutine ab2


subroutine ab2_Imp(durk1,dvrk1,dwrk1)
!
! First step of a low-storage 3rd-order Runge-Kutta scheme 
! for time integration of the momentum equations.
! Out : rhs of mom. eqs at RK1 level
!
implicit none
integer i,j,k,ii
real, intent(inout) :: durk1(0:i1,0:j1,0:k1)
real, intent(inout) :: dvrk1(0:i1,0:j1,0:k1)
real, intent(inout) :: dwrk1(0:i1,0:j1,0:k1)
real :: dnewu(0:i1,0:j1,0:k1) ! dummy array
real :: dnewv(0:i1,0:j1,0:k1) ! dummy array
real :: dneww(0:i1,0:j1,0:k1) ! dummy array
real :: duDiff(0:i1,0:j1,0:k1) ! dummy array
real :: dvDiff(0:i1,0:j1,0:k1) ! dummy array
real :: dwDiff(0:i1,0:j1,0:k1) ! dummy array
real :: durkad(0:i1,0:j1,0:k1) ! dummy array
real :: dvrkad(0:i1,0:j1,0:k1) ! dummy array
real :: dwrkad(0:i1,0:j1,0:k1) ! dummy array
real :: lapU(0:i1,0:j1,0:k1) ! dummy array
real :: lapV(0:i1,0:j1,0:k1) ! dummy array
real :: lapW(0:i1,0:j1,0:k1) ! dummy array
real :: factor,alpha,gama,zeta
real ::dPFM_couple_xnew(0:i1,0:j1,0:k1)
real ::dPFM_couple_ynew(0:i1,0:j1,0:k1)
real ::dPFM_couple_znew(0:i1,0:j1,0:k1)
real::factor1,factor2
real ,dimension(0:i1,0:j1,0:k1) ::u_ast,v_ast,w_ast
real::nu_m
nu_m =0.5*(vis_1/rho_1+vis_2/rho_2)
u_ast    = uo !2.*uo-uoo
#ifdef TwoD
u_ast=0.
#endif
v_ast    = vo !2.*vo-voo
w_ast    = wo !2.*wo-woo

factor1 = dt*(1.5)
factor2 = dt*(-0.5)
durkad                   =0.
dvrkad                   =0.
dwrkad                   =0.
duDiff                   =0.
dvDiff                   =0.
dwDiff                   =0.
dPFM_couple_xnew         =0.
dPFM_couple_ynew         =0.
dPFM_couple_znew         =0.

call momad(durkad,dvrkad,dwrkad,u_ast,v_ast,w_ast)
call momdiff(duDiff,dvDiff,dwDiff,u_ast,v_ast,w_ast)
call PFM_couple(u_ast,v_ast,w_ast,dPFM_couple_xnew,dPFM_couple_ynew,dPFM_couple_znew)
call LapVel(lapU,lapV,lapW,u_ast,v_ast,w_ast)
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      dnewu(i,j,k) = durkad(i,j,k)+dPFM_couple_xnew(i,j,k)
      dnewv(i,j,k) = dvrkad(i,j,k)+dPFM_couple_ynew(i,j,k)
      dneww(i,j,k) = dwrkad(i,j,k)+dPFM_couple_znew(i,j,k)
      dudt(i,j,k)  = unew(i,j,k) + factor1*dnewu(i,j,k) + factor2*durk1(i,j,k)+dt*duDiff(i,j,k)-0.5*nu_m*dt*lapU(i,j,k)
      dvdt(i,j,k)  = vnew(i,j,k) + factor1*dnewv(i,j,k) + factor2*dvrk1(i,j,k)+dt*dvDiff(i,j,k)-0.5*nu_m*dt*lapV(i,j,k)
      dwdt(i,j,k)  = wnew(i,j,k) + factor1*dneww(i,j,k) + factor2*dwrk1(i,j,k)+dt*dwDiff(i,j,k)-0.5*nu_m*dt*lapW(i,j,k)
      durk1(i,j,k) = dnewu(i,j,k)
      dvrk1(i,j,k) = dnewv(i,j,k)
      dwrk1(i,j,k) = dneww(i,j,k)
    enddo
  enddo
enddo


#ifdef TwoD
dudt(:,:,:) = 0.
unew(:,:,:) = 0.
#endif
return
end subroutine ab2_Imp


subroutine updatePfmDyn
implicit none
integer i,j,k,ii
real :: factor1,factor2
real ::dPFMnew(0:i1,0:j1,0:k1)
real:: Phi_wall_n(0:i1,0:j1),phi_wall_n_top(0:i1,0:j1)
real::phi_wall_nIBM(0:i1,0:j1,0:k1)
real:: dPFM_bound_new,Phi_wall_n1,Phi_wall_n1_top
real, dimension(0:i1,0:j1) ::dPFM_bound_n,dPFM_bound_n_top
real:: nx,ny,nz,n_abs,phi_m
integer:: i_p, j_p, k_p 
real:: eps = 1e-12
integer:: i_mir,j_mir,k_mir
real:: xxx,yyy,zzz
real::dn_1,dn_2,alpha
factor1 = dt*( 1.5)
factor2 = dt*(-0.5)
   do k=1,kmax
    do j=1,jmax
     do i=1,imax
         phi_wall_nIBM(i,j,k) = 0.
         nx =  nx_surf(i,j,k)
         ny =  ny_surf(i,j,k)
         nz =  nz_surf(i,j,k)
         n_abs = nabs_surf(i,j,k)
         dn_1 = deltan(i,j,k)
         dn_2 = march_step
         alpha= dn_2/(dn_1+1e-12)

         if  ( abs(nabs_surf(i,j,k)).gt.1e-12) then
               call interpolation_mirror(PFM_phi,i,j,k,phi_m)
               xxx=(i+coords(1)*imax)*dx-0.5*dx
               yyy=(j+coords(2)*jmax)*dy-0.5*dy
               zzz=(k               )*dz-0.5*dz
               phi_wall_nIBM(i,j,k) = (alpha*PFM_phi(i,j,k)+phi_m)/(alpha+1)
         endif
     enddo
    enddo
   enddo
   !update the boundary to n+1
   do k=1,kmax
    do j=1,jmax
     do i=1,imax
         ! Inside the wall
         if  ( abs(nabs_surf(i,j,k)).gt.1e-12) then
           nx =  nx_surf(i,j,k)
           ny =  ny_surf(i,j,k)
           nz =  nz_surf(i,j,k)
           n_abs = sqrt(nx*nx+ny*ny+nz*nz)
           Phi_wall_n1=phi_wall_nIBM(i,j,k) +factor1*dPFM_boundIBM(i,j,k)+factor2*dPFM_boundIBM_old(i,j,k)

           if (abs(PFM_mu_f).gt.1e-12) then
              xxx=(i+coords(1)*imax)*dx-0.5*dx
              yyy=(j+coords(2)*jmax)*dy-0.5*dy
              zzz=(k               )*dz-0.5*dz
              call interpolation_mirror(PFM_phi,i,j,k,phi_m)
              PFM_phi(i,j,k) = ((alpha+1.)*Phi_wall_n1-phi_m)/(alpha)
              dPFM_boundIBM_old(i,j,k)=dPFM_boundIBM(i,j,k)
           endif
         endif
     enddo
    enddo
   enddo


end subroutine updatePfmDyn
end module mod_rk
