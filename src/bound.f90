module mod_bound
use mod_common_mpi
implicit none
private

public bounduvw, boundp, updthalos, updthalosBig,updthalosBigInt,boundPFM,boundChem,boundloadd,boundloadIBM
contains
subroutine bounduvw(UU,VV,WW,PFM_phiB,prediction)
use mod_param
use mod_common_IBM
use mod_IBM
implicit none
integer :: i,j,k
real, dimension(-5:i1+5,-5:j1+5,-5:k1+5), intent(in):: PFM_phiB
integer,intent(in):: prediction
real:: Phi_grad_abs
real, dimension(0:i1,0:j1,0:k1), intent(inout) :: UU,VV,WW
real, dimension(0:i1,0:j1,0:k1):: Vt_ghost,Vn_ghost
real, dimension(-5:i1+5,-5:j1+5,-5:k1+5):: Uin,Vin,Win,Ucent,Vcent,Wcent
real ::  second_term_a,second_term_b,second_term
real :: mu_bound
real :: phi_wall , phi_m
real:: phi_jp, phi_jm,phi_kp,phi_km
real:: g_prime_phi,v_wall,w_wall
real:: xxx,yyy,zzz
real:: eps = 1e-12
real:: dn,x_mir,y_mir,z_mir,U_m,V_m,W_m
real:: DphiDy,DphiDz,DphiDn,DphiDx_t,v_tan
real:: DphiDy_w,DphiDz_w,Vt_m,Vn_m,Vt_g
real:: DphiDy_g,DphiDz_g,Dphidx_m,DphiDy_m,DphiDz_m
real :: Vtan_at_wall,Vtan_g,dVtdn
real:: alpha,dn_1,dn_2,numer,denum
real:: ny , nz
#ifdef NIBM
  do j=0,j1
    do i=0,i1
      UU(i,j,0)    = -UU(i,j,1)     
      WW(i,j,0)    =  0.0 
      UU(i,j,k1)   = -UU(i,j,kmax)  
      WW(i,j,kmax) =  0.0           
      WW(i,j,k1)   =  WW(i,j,kmax-1)  
      select case (PhaseField)
        case('PhasSep')
        VV(i,j,0)   = -VV(i,j,1)
        VV(i,j,k1)   = -VV(i,j,kmax)     
        case ('Droplet')
        VV(i,j,k1)   = -VV(i,j,kmax)
#ifdef SI
        VV(i,j,0)    = -VV(i,j,1)
#endif
     end select
    enddo
  enddo

    call updthalos(UU,1)
    call updthalos(VV,1)
    call updthalos(WW,1)

    call updthalos(UU,2)
    call updthalos(VV,2)
    call updthalos(WW,2)
#ifdef EXP
!Top Wall
   select case (PhaseField)
    case ('Couette')
      do j=0,j1
        do i=0,i1
          phi_jp = 0.5*(PFM_phiB(i,j+1,k1)+PFM_phiB(i,j+1,kmax))
          phi_jm = 0.5*(PFM_phiB(i,j,k1)+PFM_phiB(i,j,kmax))

          phi_kp = 0.5*(PFM_phiB(i,j,kmax)+PFM_phiB(i,j+1,kmax))
          phi_km = 0.5*(PFM_phiB(i,j,k1)+PFM_phiB(i,j+1,k1))

          phi_wall = 0.25*(PFM_phiB(i,j,k1)+PFM_phiB(i,j+1,k1)+PFM_phiB(i,j+1,kmax)+PFM_phiB(i,j,kmax))
          g_prime_phi = 0.75*(1-phi_wall*phi_wall)



          mu_bound= 0.5*(phi_wall+1)*vis_2-0.5*(phi_wall-1)*vis_1

          DphiDz = (phi_kp-phi_km)*dzi
          DphiDy = (phi_jp-phi_jm)*dyi

          second_term_a = -PFM_sigma_ref*(PFM_l)*DphiDz* 3./(2.*sqrt(2.))
          second_term_b = -PFM_sigma_ref*cos(PFM_thetta)*g_prime_phi
          second_term = (second_term_a+second_term_b)*(DphiDy)
          VV(i,j,k1)    = ((slip_length*dzi-0.5)*VV(i,j,kmax)+(second_term*slip_length/mu_bound)+vel_wall_top)/(0.5+slip_length*dzi)

        enddo
      enddo
     end select
!   !Bottom wall
    do j=0,j1
      do i=0,i1
        phi_jp = 0.5*(PFM_phiB(i,j+1,0)+PFM_phiB(i,j+1,1))
        phi_jm = 0.5*(PFM_phiB(i,j,0)+PFM_phiB(i,j,1))

        phi_kp = 0.5*(PFM_phiB(i,j,1)+PFM_phiB(i,j+1,1))
        phi_km = 0.5*(PFM_phiB(i,j,0)+PFM_phiB(i,j+1,0))

        phi_wall = 0.25*(PFM_phiB(i,j,0)+PFM_phiB(i,j+1,0)+PFM_phiB(i,j+1,1)+PFM_phiB(i,j,1))
        g_prime_phi = 0.75*(1-phi_wall*phi_wall) 

        mu_bound= 0.5*(phi_wall+1)*vis_2-0.5*(phi_wall-1)*vis_1

        DphiDz = (phi_kp-phi_km)*dzi
        DphiDy = (phi_jp-phi_jm)*dyi

        second_term_a = -PFM_sigma_ref*(PFM_l)*DphiDz* 3./(2.*sqrt(2.))
        second_term_b = -PFM_sigma_ref*cos(PFM_thetta)*g_prime_phi
        second_term = (second_term_a+second_term_b)*(DphiDy)
        VV(i,j,0)    = ((slip_length*dzi-0.5)*VV(i,j,1)+(second_term*slip_length/mu_bound)+vel_wall_bot)/(0.5+slip_length*dzi)

      enddo
    enddo

! communicate data in x direction (periodic b.c.'s incorporated)
call updthalos(UU,1)
call updthalos(VV,1)
call updthalos(WW,1)
! communicate data in y direction (periodic b.c.'s incorporated)
call updthalos(UU,2)
call updthalos(VV,2)
call updthalos(WW,2)
#endif
#endif
#ifdef IBM
 do j=0,j1
   do i=0,i1
      UU(i,j,k1)   = -UU(i,j,kmax)  
      VV(i,j,k1)   = -VV(i,j,kmax)    
      WW(i,j,kmax) = 0.             
      WW(i,j,k1)   = WW(i,j,kmax-1)  

      UU(i,j,0)   = -UU(i,j,1)  
      VV(i,j,0)   = -VV(i,j,1)    
      WW(i,j,0) = 0.              
    enddo
  enddo



    call updthalos(UU,1)
    call updthalos(VV,1)
    call updthalos(WW,1)

    call updthalos(UU,2)
    call updthalos(VV,2)
    call updthalos(WW,2)


    do  k = 0,k1
      do j = 0,j1
       do i = 0,i1
        Uin(i,j,k) = UU(i,j,k)
        Vin(i,j,k) = VV(i,j,k)
        Win(i,j,k) = WW(i,j,k)
  
        Uin(i,j,-1  ) = UU(i,j,0)
        Uin(i,j,-2  ) = UU(i,j,0)
        Uin(i,j,-3  ) = UU(i,j,0)
        Uin(i,j,-4  ) = UU(i,j,0)
        Uin(i,j,-5  ) = UU(i,j,0)
        Uin(i,j,k1+1) = UU(i,j,k1)
        Uin(i,j,k1+2) = UU(i,j,k1)
        Uin(i,j,k1+3) = UU(i,j,k1)
        Uin(i,j,k1+4) = UU(i,j,k1)
        Uin(i,j,k1+5) = UU(i,j,k1)
  
        Vin(i,j,-1)   = VV(i,j,0)
        Vin(i,j,-2)   = VV(i,j,0)
        Vin(i,j,-3)   = VV(i,j,0)
        Vin(i,j,-4)   = VV(i,j,0)
        Vin(i,j,-5)   = VV(i,j,0)
        Vin(i,j,k1+1) = VV(i,j,k1)
        Vin(i,j,k1+2) = VV(i,j,k1)
        Vin(i,j,k1+3) = VV(i,j,k1)
        Vin(i,j,k1+4) = VV(i,j,k1)
        Vin(i,j,k1+5) = VV(i,j,k1)
  
        Win(i,j,-1)   = WW(i,j,0)
        Win(i,j,-2)   = WW(i,j,0)
        Win(i,j,-3)   = WW(i,j,0)
        Win(i,j,-4)   = WW(i,j,0)
        Win(i,j,-5)   = WW(i,j,0)
        Win(i,j,k1+1) = WW(i,j,k1)
        Win(i,j,k1+2) = WW(i,j,k1)
        Win(i,j,k1+3) = WW(i,j,k1)
        Win(i,j,k1+4) = WW(i,j,k1)
        Win(i,j,k1+5) = WW(i,j,k1)
       enddo
      enddo
    enddo

    call updthalosbig(Uin,1)
    call updthalosbig(Uin,2)
    call updthalosbig(Vin,1)
    call updthalosbig(Vin,2)
    call updthalosbig(Win,1)
    call updthalosbig(Win,2)

! Cell center Velocity
    do  k = 0,k1
      do j = 0,j1
       do i = 0,i1
        Ucent(i,j,k) = 0.5*(Uin(i,j,k)+Uin(i-1,j,k))
        Vcent(i,j,k) = 0.5*(Vin(i,j,k)+Vin(i,j-1,k))
        Wcent(i,j,k) = 0.5*(Win(i,j,k)+Win(i,j,k-1))
     
        Ucent(i,j,-1  ) = 0.5*(Uin(i,j,0)+Uin(i-1,j,0)) 
        Ucent(i,j,-2  ) = 0.5*(Uin(i,j,0)+Uin(i-1,j,0))
        Ucent(i,j,-3  ) = 0.5*(Uin(i,j,0)+Uin(i-1,j,0))
        Ucent(i,j,-4  ) = 0.5*(Uin(i,j,0)+Uin(i-1,j,0))
        Ucent(i,j,-5  ) = 0.5*(Uin(i,j,0)+Uin(i-1,j,0))
        Ucent(i,j,k1+1) = 0.5*(Uin(i,j,k1)+Uin(i-1,j,k1))
        Ucent(i,j,k1+2) = 0.5*(Uin(i,j,k1)+Uin(i-1,j,k1))
        Ucent(i,j,k1+3) = 0.5*(Uin(i,j,k1)+Uin(i-1,j,k1))
        Ucent(i,j,k1+4) = 0.5*(Uin(i,j,k1)+Uin(i-1,j,k1))
        Ucent(i,j,k1+5) = 0.5*(Uin(i,j,k1)+Uin(i-1,j,k1))
     
        Vcent(i,j,-1)   = 0.5*(Vin(i,j,0)+Vin(i,j-1,0))
        Vcent(i,j,-2)   = 0.5*(Vin(i,j,0)+Vin(i,j-1,0))
        Vcent(i,j,-3)   = 0.5*(Vin(i,j,0)+Vin(i,j-1,0))
        Vcent(i,j,-4)   = 0.5*(Vin(i,j,0)+Vin(i,j-1,0))
        Vcent(i,j,-5)   = 0.5*(Vin(i,j,0)+Vin(i,j-1,0))
        Vcent(i,j,k1+1) = 0.5*(Vin(i,j,k1)+Vin(i,j-1,k1))
        Vcent(i,j,k1+2) = 0.5*(Vin(i,j,k1)+Vin(i,j-1,k1))
        Vcent(i,j,k1+3) = 0.5*(Vin(i,j,k1)+Vin(i,j-1,k1))
        Vcent(i,j,k1+4) = 0.5*(Vin(i,j,k1)+Vin(i,j-1,k1))
        Vcent(i,j,k1+5) = 0.5*(Vin(i,j,k1)+Vin(i,j-1,k1))
  
        Wcent(i,j,-1)   = 0.5*(Win(i,j,0)+Win(i,j,-1))
        Wcent(i,j,-2)   = 0.5*(Win(i,j,0)+Win(i,j,-1))
        Wcent(i,j,-3)   = 0.5*(Win(i,j,0)+Win(i,j,-1))
        Wcent(i,j,-4)   = 0.5*(Win(i,j,0)+Win(i,j,-1))
        Wcent(i,j,-5)   = 0.5*(Win(i,j,0)+Win(i,j,-1))
        Wcent(i,j,k1+1) = 0.5*(Win(i,j,kmax)+Win(i,j,k1))
        Wcent(i,j,k1+2) = 0.5*(Win(i,j,kmax)+Win(i,j,k1))
        Wcent(i,j,k1+3) = 0.5*(Win(i,j,kmax)+Win(i,j,k1))
        Wcent(i,j,k1+4) = 0.5*(Win(i,j,kmax)+Win(i,j,k1))
        Wcent(i,j,k1+5) = 0.5*(Win(i,j,kmax)+Win(i,j,k1))
       enddo
      enddo
    enddo
    call updthalosbig(Ucent,1)
    call updthalosbig(Ucent,2)
    call updthalosbig(Vcent,1)
    call updthalosbig(Vcent,2)
    call updthalosbig(Wcent,1)
    call updthalosbig(Wcent,2)
    do k = 1, kmax
      do j=1,jmax
        do i=1,imax
         vt_ghost(i,j,k)     = 0.
         Vt_g                = 0
         if  ( abs(nabs_surf(i,j,k)).gt.1e-12) then
             ! Initializations
              DphiDn=0.
              v_wall=0.
              w_wall=0.
              second_term=0.
              xxx=(i+coords(1)*imax)*dx-0.5*dx
              yyy=(j+coords(2)*jmax)*dy-0.5*dy
              zzz=(k               )*dz-0.5*dz

              ny = ny_surf(i,j,k)
              nz = nz_surf(i,j,k)

             ! Geometrical requirements
              dn_1 = deltan(i,j,k) 
              dn_2 = march_step 
              alpha= dn_2/(dn_1+1e-12)
              if (dn_1.lt.1e-11) alpha = 1e-9
              if (dn_1.gt.1e11) alpha = 1e9              
             ! Calculating the tangential velocity at the wall
              call interpolation_2D_velocity(Ucent,Vcent,Wcent,i,j,k,U_m,V_m,W_m)
              v_wall =  (alpha*Vcent(i,j,k)+V_m)/(alpha+1)
              w_wall =  (alpha*Wcent(i,j,k)+W_m)/(alpha+1)
              v_tan =  nz*v_wall-ny*w_wall
             ! Calculating DphiDn at the wall
              call interpolation_mirror(PFM_phiB,i,j,k,phi_m)
              phi_wall = (alpha*PFM_phiB(i,j,k)+phi_m)/(alpha+1)
              mu_bound=0.5*(phi_wall+1)*vis_2-0.5*(phi_wall-1)*vis_1
              g_prime_phi = 0.75* (1-phi_wall*phi_wall)

              numer = phi_m+(alpha*alpha-1)*phi_wall-alpha*alpha*PFM_phiB(i,j,k)
              denum = alpha*(alpha+1)*dn_1
              DphiDn=numer/denum 
              ! Interpolating DphiDy and DphiDz at the wall
              call interpolation_2D_dphi(i,j,k,DphiDx_m,DphiDy_m,DphiDz_m)
              phi_jp = 0.5*(PFM_phiB(i,j,k)+PFM_phiB(i,j+1,k))
              phi_jm = 0.5*(PFM_phiB(i,j,k)+PFM_phiB(i,j-1,k))
              phi_kp = 0.5*(PFM_phiB(i,j,k)+PFM_phiB(i,j,k+1))
              phi_km = 0.5*(PFM_phiB(i,j,k)+PFM_phiB(i,j,k-1))
              DphiDy_g = (phi_jp-phi_jm)*dyi
              DphiDz_g = (phi_kp-phi_km)*dzi
              DphiDy_w = (alpha*DphiDy_g+DphiDy_m)/(alpha+1)
              DphiDz_w = (alpha*DphiDz_g+DphiDz_m)/(alpha+1)
             ! Calculating DphiDxt (tangential component)
              DphiDx_t = nz*DphiDy_w-ny*DphiDz_w
              Vt_m     = nz*V_m - ny*W_m
              Vn_m     = ny*V_m + nz*W_m
              Vtan_g   = nz*Vcent(i,j,k)-ny*Wcent(i,j,k)
              dVtdn    = (Vt_m+(alpha*alpha-1)*v_tan-alpha*alpha*Vtan_g)/(alpha*(alpha+1)*dn_1)
           ! Computing slip velocity
              second_term_a= (-PFM_sigma_ref*PFM_l)* (DphiDn)*3./(2.*sqrt(2.))
              second_term_b= (-PFM_sigma_ref)*cos(PFM_thetta)*g_prime_phi
              second_term = (second_term_a +second_term_b)*DphiDx_t!*(abs(PFM_mu_f)/(abs(PFM_mu_f)+1e-12))
            ! Calculating new tangential velocity at ghost point
              numer = Vt_m*(slip_length-dn_1)+(slip_length/mu_bound)*(alpha+1)*dn_1*second_term
              denum = alpha*dn_1+slip_length
              Vt_g  =  numer /denum
              Vt_ghost(i,j,k) = Vt_g
              Vn_ghost(i,j,k) = -Vn_m/alpha

        endif
      enddo
     enddo
   enddo
   
   
   do k = 1, kmax
      do j=1,jmax
        do i=1,imax
         ! Ghost points  :
          ny = ny_surf(i,j,k)
          nz = nz_surf(i,j,k)
          if  ( abs(nabs_surf(i,j,k)).gt.1e-12) then
           Vcent(i,j,k) =  nz*Vt_ghost(i,j,k)+ ny*Vn_ghost(i,j,k)
           Wcent(i,j,k) = -ny*Vt_ghost(i,j,k)+nz*Vn_ghost(i,j,k)
          endif
        enddo
     enddo
   enddo

if (slip_length.gt.1.e-12) then
   if (prediction.eq.1) then
     call Penalization_center(Ucent,Vcent,Wcent)
   endif ! !prediction = 

   call updthalosbig(Ucent,1)
   call updthalosbig(Ucent,2)
   call updthalosbig(Vcent,1)
   call updthalosbig(Vcent,2)
   call updthalosbig(Wcent,1)
   call updthalosbig(Wcent,2)


    do k=0,k1
      do j=0,j1
        do i=0,i1
         if( (cell_u_tag(i,j,k).lt.1).and.(cell_u_tag(i,j,k).gt.1e-12)) then
         UU(i,j,k) = 0.5*(Ucent(i,j,k)+Ucent(i+1,j,k))
         endif
         if( (cell_v_tag(i,j,k).lt.1).and.(cell_v_tag(i,j,k).gt.1e-12)) then
         VV(i,j,k) = 0.5*(Vcent(i,j,k)+Vcent(i,j+1,k))
         endif
         if( (cell_w_tag(i,j,k).lt.1).and.(cell_w_tag(i,j,k).gt.1e-12)) then
         WW(i,j,k) = 0.5*(Wcent(i,j,k)+Wcent(i,j,k+1))
         endif
         enddo
      enddo
    enddo
endif 
    if (prediction.eq.1) then
         call Penalization_face(UU,VV,WW)
    endif
    
       call updthalos(UU,1)
       call updthalos(UU,2)
       call updthalos(VV,1)
       call updthalos(VV,2)
       call updthalos(WW,1)
       call updthalos(WW,2)
#endif
return
end subroutine bounduvw
!
subroutine boundp(p)
use mod_param
use mod_common_IBM

implicit none
integer :: i,j
real, dimension(0:i1,0:j1,0:k1),intent(inout) :: p
do j=0,j1
  do i=0,i1
    p(i,j,0) = p(i,j,1)     ! Newmann (consistent with no/free-slip)
    p(i,j,k1) = p(i,j,kmax) ! Newmann (consistent with no/free-slip)
  enddo
enddo


call updthalos(p,1)
!
call updthalos(p,2)
!
return
end subroutine boundp
!******************* Phase Field Method ********************************************************
subroutine boundPFM(PFM_phiB,dPFM_boundd,dPFM_boundd_top,u,v,w)
use mod_param
use mod_interface
use mod_common_IBM
use mod_IBM
implicit none
integer :: i,j,k
real, dimension(-5:i1+5,-5:j1+5,-5:k1+5), intent(inout) :: PFM_phiB
real, dimension(0:i1,0:j1,0:k1), intent(in) :: u,v,w
real, dimension(0:i1,0:j1), intent(out) ::dPFM_boundd,dPFM_boundd_top
real:: first_term,second_term,third_term,DphiDy,DphiDz
real:: g_prime_phi,v_wall
real::phi_jp,phi_jm,phi_wall
real::PFM_mu_f_1
real, dimension(-5:i1+5,-5:j1+5,-5:k1+5) :: Uin,Vin,Win,Ucent,Vcent,Wcent
real:: dn,U_m,V_m,W_m
real:: xxx,yyy,zzz
real:: DphiDn
real::phi_m,w_wall
real::phi_kp,phi_km
real::DphiDy_g,DphiDz_g,DphiDy_w,DphiDz_w,DphiDx_m,DphiDy_m,DphiDz_m,DphiDx_tan
real::v_tan,eps = 1e-12
real:: alpha,dn_1,dn_2,numer,denum
real:: ny, nz
#ifdef NIBM
#ifdef EXP 
 select case (PhaseField)
    case('PhasSep')
      do j=-5,j1+5
        do i=-5,i1+5
          PFM_phiB(i,j,-4:0) = PFM_phiB(i,j,1)
          PFM_phiB(i,j,k1:k1+4) = PFM_phiB(i,j,kmax)
        enddo
      enddo
    case ('Droplet')
      do j=-5,j1+5
        do i=-5,i1+5
          PFM_phiB(i,j,k1:k1+4) = PFM_phiB(i,j,kmax) 
        enddo
      enddo
      dPFM_boundd_top(:,:) = 0.
    case ('Couette')
   
   !Top Wall
      do j=1,jmax
        do i=1,imax
         phi_wall = 0.5*(PFM_phiB(i,j,k1)+PFM_phiB(i,j,kmax))
         g_prime_phi = 0.75* (1-phi_wall*phi_wall)
   
         if (abs(PFM_mu_f).gt.1e-12) then
            PFM_mu_f_1 =PFM_mu_f
            dPFM_boundd_top(i,j)= 0.0
            DphiDy=0.
            DphiDz=0.
            v_wall=0.
            first_term= 0.
            second_term=0.
            third_term=0.
   
            v_wall = 0.25*( v(i,j,kmax)+v(i,j,k1)+v(i,j-1,kmax)+v(i,j-1,k1) )
            phi_jp = 0.25*(PFM_phiB(i,j,k1)+PFM_phiB(i,j+1,k1)+PFM_phiB(i,j+1,kmax)+PFM_phiB(i,j,kmax))
            phi_jm = 0.25*(PFM_phiB(i,j,k1)+PFM_phiB(i,j-1,k1)+PFM_phiB(i,j-1,kmax)+PFM_phiB(i,j,kmax))
   
            DphiDy = (phi_jp-phi_jm)*dyi
            DphiDz =(PFM_phiB(i,j,kmax)-PFM_phiB(i,j,k1))*dzi
   
            first_term = -v_wall*DphiDy
            second_term=(PFM_sigma_ref/PFM_mu_f_1)* DphiDz* 3./(2.*sqrt(2.))
            third_term= (PFM_sigma_ref/(PFM_l*PFM_mu_f_1))*cos(PFM_thetta)*g_prime_phi
            dPFM_boundd_top(i,j)= first_term+second_term+third_term
         else
            PFM_phiB(i,j,k1) = PFM_phiB(i,j,kmax)+ (dz/PFM_l)*g_prime_phi*cos(PFM_thetta)*(2.*sqrt(2.)/3)
            dPFM_boundd_top(i,j) = 0.

          endif
        enddo
      enddo
   
   end select
   !Bottom wall
      do j=1,jmax
       do i=1,imax
         phi_wall = 0.5*(PFM_phiB(i,j,0)+PFM_phiB(i,j,1))
         g_prime_phi = 0.75* (1-phi_wall*phi_wall)
         
         if (abs(PFM_mu_f).gt.1e-12) then
            PFM_mu_f_1 =PFM_mu_f
            dPFM_boundd(i,j)= 0.0
            DphiDy=0.
            DphiDz=0.
            v_wall=0.
            first_term= 0.
            second_term=0.
            third_term=0.
   
            v_wall = 0.25*( v(i,j,1)+v(i,j,0)+v(i,j-1,1)+v(i,j-1,0) )
            phi_jp = 0.25*(PFM_phiB(i,j,0)+PFM_phiB(i,j+1,0)+PFM_phiB(i,j+1,1)+PFM_phiB(i,j,1))
            phi_jm = 0.25*(PFM_phiB(i,j,0)+PFM_phiB(i,j-1,0)+PFM_phiB(i,j-1,1)+PFM_phiB(i,j,1))
   
            DphiDy = (phi_jp-phi_jm)*dyi
            DphiDz =(PFM_phiB(i,j,1)-PFM_phiB(i,j,0))*dzi
   
            first_term = -v_wall*DphiDy
            second_term=(PFM_sigma_ref/PFM_mu_f_1)* DphiDz* 3./(2.*sqrt(2.))
            third_term= (PFM_sigma_ref/(PFM_l*PFM_mu_f_1))*cos(PFM_thetta)*g_prime_phi
            dPFM_boundd(i,j)= first_term+second_term+third_term
         else
              PFM_phiB(i,j,0) = PFM_phiB(i,j,1)+ (dz/PFM_l)*g_prime_phi*cos(PFM_thetta)*(2.*sqrt(2.)/3)
              dPFM_boundd(i,j)= 0.  
        endif 
        enddo
      enddo
#endif
#endif
#ifdef IBM
  do k=0,k1
   do j=0,j1
     do i=0,i1
       Uin(i,j,k) = u(i,j,k)
       Vin(i,j,k) = v(i,j,k)
       Win(i,j,k) = w(i,j,k)
       Uin(i,j,k) = u(i,j,k)
       Vin(i,j,k) = v(i,j,k)
       Win(i,j,k) = w(i,j,k)
    
       Uin(i,j,-1  ) = u(i,j,0)
       Uin(i,j,-2  ) = u(i,j,0)
       Uin(i,j,-3  ) = u(i,j,0)
       Uin(i,j,-4  ) = u(i,j,0)
       Uin(i,j,-5  ) = u(i,j,0)
       Uin(i,j,k1+1) = u(i,j,k1)
       Uin(i,j,k1+2) = u(i,j,k1)
       Uin(i,j,k1+3) = u(i,j,k1)
       Uin(i,j,k1+4) = u(i,j,k1)
       Uin(i,j,k1+5) = u(i,j,k1)
    
       Vin(i,j,-1)   = v(i,j,0)
       Vin(i,j,-2)   = v(i,j,0)
       Vin(i,j,-3)   = v(i,j,0)
       Vin(i,j,-4)   = v(i,j,0)
       Vin(i,j,-5)   = v(i,j,0)
       Vin(i,j,k1+1) = v(i,j,k1)
       Vin(i,j,k1+2) = v(i,j,k1)
       Vin(i,j,k1+3) = v(i,j,k1)
       Vin(i,j,k1+4) = v(i,j,k1)
       Vin(i,j,k1+5) = v(i,j,k1)
    
       Win(i,j,-1)   = w(i,j,0)
       Win(i,j,-2)   = w(i,j,0)
       Win(i,j,-3)   = w(i,j,0)
       Win(i,j,-4)   = w(i,j,0)
       Win(i,j,-5)   = w(i,j,0)
       Win(i,j,k1+1) = w(i,j,k1)
       Win(i,j,k1+2) = w(i,j,k1)
       Win(i,j,k1+3) = w(i,j,k1)
       Win(i,j,k1+4) = w(i,j,k1)
       Win(i,j,k1+5) = w(i,j,k1)
      enddo
     enddo
   enddo
   call updthalosBig(Uin,1)
   call updthalosBig(Uin,2)
   call updthalosBig(Vin,1)
   call updthalosBig(Vin,2)
   call updthalosBig(Win,1)
   call updthalosBig(Win,2)

! Cell center Velocity
    do  k = 0,k1
      do j = 0,j1
       do i = 0,i1
        Ucent(i,j,k) = 0.5*(Uin(i,j,k)+Uin(i-1,j,k))
        Vcent(i,j,k) = 0.5*(Vin(i,j,k)+Vin(i,j-1,k))
        Wcent(i,j,k) = 0.5*(Win(i,j,k)+Win(i,j,k-1))
     
        Ucent(i,j,-1  ) = 0.5*(Uin(i,j,0)+Uin(i-1,j,0)) 
        Ucent(i,j,-2  ) = 0.5*(Uin(i,j,0)+Uin(i-1,j,0))
        Ucent(i,j,-3  ) = 0.5*(Uin(i,j,0)+Uin(i-1,j,0))
        Ucent(i,j,-4  ) = 0.5*(Uin(i,j,0)+Uin(i-1,j,0))
        Ucent(i,j,-5  ) = 0.5*(Uin(i,j,0)+Uin(i-1,j,0))
        Ucent(i,j,k1+1) = 0.5*(Uin(i,j,k1)+Uin(i-1,j,k1))
        Ucent(i,j,k1+2) = 0.5*(Uin(i,j,k1)+Uin(i-1,j,k1))
        Ucent(i,j,k1+3) = 0.5*(Uin(i,j,k1)+Uin(i-1,j,k1))
        Ucent(i,j,k1+4) = 0.5*(Uin(i,j,k1)+Uin(i-1,j,k1))
        Ucent(i,j,k1+5) = 0.5*(Uin(i,j,k1)+Uin(i-1,j,k1))
     
        Vcent(i,j,-1)   = 0.5*(Vin(i,j,0)+Vin(i,j-1,0))
        Vcent(i,j,-2)   = 0.5*(Vin(i,j,0)+Vin(i,j-1,0))
        Vcent(i,j,-3)   = 0.5*(Vin(i,j,0)+Vin(i,j-1,0))
        Vcent(i,j,-4)   = 0.5*(Vin(i,j,0)+Vin(i,j-1,0))
        Vcent(i,j,-5)   = 0.5*(Vin(i,j,0)+Vin(i,j-1,0))
        Vcent(i,j,k1+1) = 0.5*(Vin(i,j,k1)+Vin(i,j-1,k1))
        Vcent(i,j,k1+2) = 0.5*(Vin(i,j,k1)+Vin(i,j-1,k1))
        Vcent(i,j,k1+3) = 0.5*(Vin(i,j,k1)+Vin(i,j-1,k1))
        Vcent(i,j,k1+4) = 0.5*(Vin(i,j,k1)+Vin(i,j-1,k1))
        Vcent(i,j,k1+5) = 0.5*(Vin(i,j,k1)+Vin(i,j-1,k1))
  
        Wcent(i,j,-1)   = 0.5*(Win(i,j,0)+Win(i,j,-1))
        Wcent(i,j,-2)   = 0.5*(Win(i,j,0)+Win(i,j,-1))
        Wcent(i,j,-3)   = 0.5*(Win(i,j,0)+Win(i,j,-1))
        Wcent(i,j,-4)   = 0.5*(Win(i,j,0)+Win(i,j,-1))
        Wcent(i,j,-5)   = 0.5*(Win(i,j,0)+Win(i,j,-1))
        Wcent(i,j,k1+1) = 0.5*(Win(i,j,kmax)+Win(i,j,k1))
        Wcent(i,j,k1+2) = 0.5*(Win(i,j,kmax)+Win(i,j,k1))
        Wcent(i,j,k1+3) = 0.5*(Win(i,j,kmax)+Win(i,j,k1))
        Wcent(i,j,k1+4) = 0.5*(Win(i,j,kmax)+Win(i,j,k1))
        Wcent(i,j,k1+5) = 0.5*(Win(i,j,kmax)+Win(i,j,k1))
       enddo
      enddo
    enddo
    call updthalosbig(Ucent,1)
    call updthalosbig(Ucent,2)
    call updthalosbig(Vcent,1)
    call updthalosbig(Vcent,2)
    call updthalosbig(Wcent,1)
    call updthalosbig(Wcent,2)



   do k = 1, kmax
      do j=1,jmax
        do i=1,imax
          ! Ghost points  :
          if  ( abs(nabs_surf(i,j,k)).gt.1e-12) then
                call interpolation_mirror(PFM_phiB,i,j,k,phi_m)
                PFM_mu_f_1 =PFM_mu_f
                dPFM_boundIBM(i,j,k)= 0.0
                DphiDn=0.
                v_wall=0.
                w_wall=0.
                ny = ny_surf(i,j,k)
                nz = nz_surf(i,j,k)
                first_term= 0.
                second_term=0.
                third_term=0.
                xxx=(i+coords(1)*imax)*dx-0.5*dx
                yyy=(j+coords(2)*jmax)*dy-0.5*dy
                zzz=(k               )*dz-0.5*dz
                v_wall = 0.
                w_wall = 0.
                dn_1 = deltan(i,j,k)
                dn_2 = march_step 
                alpha= dn_2/(dn_1+1e-12)
                if (dn_1.lt.1e-11) alpha = 1e-9
                if (dn_1.gt.1e11) alpha = 1e9
                phi_wall = (alpha*PFM_phiB(i,j,k)+phi_m)/(alpha+1)
                g_prime_phi = 0.75* (1.0-phi_wall*phi_wall)
               !Dynamic BC
                if (abs(PFM_mu_f).gt.1e-12) then
                    numer = phi_m+(alpha*alpha-1)*phi_wall-alpha*alpha*PFM_phiB(i,j,k)
                    denum = alpha*(alpha+1)*dn_1
                    DphiDn =numer/denum                    
                    call interpolation_2D_velocity(Ucent,Vcent,Wcent,i,j,k,U_m,V_m,W_m)
                    v_wall =  (alpha*Vcent(i,j,k)+V_m)/(alpha+1)
                    w_wall =  (alpha*Wcent(i,j,k)+W_m)/(alpha+1)
!*****************************************************************************
                   ! Calculating the tangential velocity at the wall
                    v_tan =  nz*v_wall-ny*w_wall
                    call interpolation_2D_dphi(i,j,k,DphiDx_m,DphiDy_m,DphiDz_m)
                    phi_jp = 0.5*(PFM_phiB(i,j,k)+PFM_phiB(i,j+1,k))
                    phi_jm = 0.5*(PFM_phiB(i,j,k)+PFM_phiB(i,j-1,k))
                    phi_kp = 0.5*(PFM_phiB(i,j,k)+PFM_phiB(i,j,k+1))
                    phi_km = 0.5*(PFM_phiB(i,j,k)+PFM_phiB(i,j,k-1))
                    DphiDy_g = (phi_jp-phi_jm)*dyi
                    DphiDz_g = (phi_kp-phi_km)*dzi
                    DphiDy_w = (alpha*DphiDy_g+DphiDy_m)/(alpha+1)
                    DphiDz_w = (alpha*DphiDz_g+DphiDz_m)/(alpha+1)
                    DphiDx_tan = nz*DphiDy_w-ny*DphiDz_w
!*****************************************************************************
                    first_term = 0.
                    if (abs(slip_length).gt.1e-12) first_term = -v_tan*DphiDx_tan
                    second_term= (PFM_sigma_ref/PFM_mu_f_1)* (DphiDn)*3./(2.*sqrt(2.))
                    third_term = (PFM_sigma_ref/(PFM_l*PFM_mu_f_1))*cos(PFM_thetta)*g_prime_phi
                    dPFM_boundIBM(i,j,k)=  first_term+second_term+third_term
                  !Static B,
                else
                    first_term = (2*sqrt(2.)/3)*g_prime_phi*cos(PFM_thetta)/PFM_l
                    PFM_phiB(i,j,k) =  phi_m+first_term*(alpha+1)*dn_1
                    dPFM_boundIBM(i,j,k)=  0.

                endif
            !Grids inside solid without any neighbor in fluid
             else 
                dPFM_boundIBM(i,j,k)= 0.
          endif

       enddo
     enddo
   enddo

#endif
do j=0,j1
  do i=0,i1
#ifdef IBM
  PFM_phiB(i,j,0) = PFM_phiB(i,j,1)
  PFM_phiB(i,j,k1) = PFM_phiB(i,j,kmax)
#endif
#ifdef SI
  PFM_phiB(i,j,0) = PFM_phiB(i,j,1)
  PFM_phiB(i,j,k1) = PFM_phiB(i,j,kmax)
#endif
    PFM_phiB(i,j,-1)    = PFM_phiB(i,j,0)
    PFM_phiB(i,j,-2)    = PFM_phiB(i,j,0)
    PFM_phiB(i,j,-3)    = PFM_phiB(i,j,0)
    PFM_phiB(i,j,-4)    = PFM_phiB(i,j,0)
    PFM_phiB(i,j,-5)    = PFM_phiB(i,j,0)
    PFM_phiB(i,j,k1+1)  = PFM_phiB(i,j,k1)
    PFM_phiB(i,j,k1+2)  = PFM_phiB(i,j,k1)
    PFM_phiB(i,j,k1+3)  = PFM_phiB(i,j,k1)
    PFM_phiB(i,j,k1+4)  = PFM_phiB(i,j,k1)
    PFM_phiB(i,j,k1+5)  = PFM_phiB(i,j,k1)
  enddo
enddo
call updthalosBig(PFM_phiB,1)
call updthalosBig(PFM_phiB,2)
!
return
end subroutine boundPFM

subroutine boundChem(chem_potB)
use mod_param
use mod_IBM
use mod_common_IBM
implicit none
real, dimension(-5:i1+5,-5:j1+5,-5:k1+5), intent(inout):: chem_potB
integer:: i,j,k
real::chem_m
do j=-5,j1+5
  do i=-5,i1+5
    chem_PotB(i,j,0)    = Chem_potB(i,j,1)
    chem_PotB(i,j,-1)   = Chem_potB(i,j,1)
    chem_PotB(i,j,-2)   = Chem_potB(i,j,1)
    chem_PotB(i,j,-3)   = Chem_potB(i,j,1)
    chem_PotB(i,j,-4)   = Chem_potB(i,j,1)
    chem_PotB(i,j,-5)   = Chem_potB(i,j,1)
    chem_PotB(i,j,k1+5) = Chem_potB(i,j,kmax)
    chem_PotB(i,j,k1+4) = Chem_potB(i,j,kmax)
    chem_PotB(i,j,k1+3) = Chem_potB(i,j,kmax)
    chem_PotB(i,j,k1+2) = Chem_potB(i,j,kmax)
    chem_PotB(i,j,k1+1) = Chem_potB(i,j,kmax)
    chem_PotB(i,j,k1)   = Chem_potB(i,j,kmax)
 enddo
enddo
#ifdef IBM
   do k = 1, kmax
      do j=1,jmax
        do i=1,imax
          ! Ghost points  :
          if  ( abs(nabs_surf(i,j,k)).gt.1e-12) then
            call interpolation_mirror(chem_PotB,i,j,k,chem_m)
            chem_PotB(i,j,k) = chem_m
          endif
        enddo
      enddo
   enddo
#endif
call updthalosBig(chem_PotB,1)
call updthalosBig(chem_PotB,2)
end subroutine boundChem

subroutine boundloadd(u,v,w,p,PFM_phiB)
use mod_param
implicit none
real, dimension(-5:i1+5,-5:j1+5,-5:k1+5), intent(inout) :: PFM_phiB
real, dimension(0:i1,0:j1,0:k1), intent(inout) :: u,v,w,p
integer::i,j
 do j=0,j1
   do i=0,i1
       PFM_phiB(i,j,-1)    = PFM_phiB(i,j,0)
       PFM_phiB(i,j,-2)    = PFM_phiB(i,j,0)
       PFM_phiB(i,j,-3)    = PFM_phiB(i,j,0)
       PFM_phiB(i,j,-4)    = PFM_phiB(i,j,0)
       PFM_phiB(i,j,-5)    = PFM_phiB(i,j,0)
       PFM_phiB(i,j,k1+1)  = PFM_phiB(i,j,k1)
       PFM_phiB(i,j,k1+2)  = PFM_phiB(i,j,k1)
       PFM_phiB(i,j,k1+3)  = PFM_phiB(i,j,k1)
       PFM_phiB(i,j,k1+4)  = PFM_phiB(i,j,k1)
       PFM_phiB(i,j,k1+5)  = PFM_phiB(i,j,k1)
    enddo
 enddo
call updthalos(u,1)
call updthalos(v,1)
call updthalos(w,1)
call updthalos(u,2)
call updthalos(v,2)
call updthalos(w,2)
call updthalos(p,1)
call updthalos(p,2)
call updthalosBig(PFM_phiB,1)
call updthalosBig(PFM_phiB,2)
end subroutine boundloadd

subroutine boundloadIBM(cell_phi_tagB,cell_u_tagB,cell_v_tagB,cell_w_tagB, &
                      nx_surfB,ny_surfB,nz_surfB,nabs_surfB)
use mod_param
implicit none
real, dimension(-5:i1+5,-5:j1+5,-5:k1+5), intent(inout) :: cell_phi_tagB,cell_u_tagB,cell_v_tagB
real, dimension(-5:i1+5,-5:j1+5,-5:k1+5), intent(inout) :: cell_w_tagB,nx_surfB,ny_surfB,nz_surfB,nabs_surfB
integer::i,j
 do j=0,j1
   do i=0,i1
       cell_phi_tagB(i,j,0)     = cell_phi_tagB(i,j,1)
       cell_phi_tagB(i,j,-1)    = cell_phi_tagB(i,j,1)
       cell_phi_tagB(i,j,-2)    = cell_phi_tagB(i,j,1)
       cell_phi_tagB(i,j,-3)    = cell_phi_tagB(i,j,1)
       cell_phi_tagB(i,j,-4)    = cell_phi_tagB(i,j,1)
       cell_phi_tagB(i,j,-5)    = cell_phi_tagB(i,j,1)

       cell_phi_tagB(i,j,k1)    = cell_phi_tagB(i,j,kmax)
       cell_phi_tagB(i,j,k1+1)  = cell_phi_tagB(i,j,kmax)
       cell_phi_tagB(i,j,k1+2)  = cell_phi_tagB(i,j,kmax)
       cell_phi_tagB(i,j,k1+3)  = cell_phi_tagB(i,j,kmax)
       cell_phi_tagB(i,j,k1+4)  = cell_phi_tagB(i,j,kmax)
       cell_phi_tagB(i,j,k1+5)  = cell_phi_tagB(i,j,kmax)



       cell_u_tagB(i,j,0)     = cell_u_tagB(i,j,1)
       cell_u_tagB(i,j,-1)    = cell_u_tagB(i,j,1)
       cell_u_tagB(i,j,-2)    = cell_u_tagB(i,j,1)
       cell_u_tagB(i,j,-3)    = cell_u_tagB(i,j,1)
       cell_u_tagB(i,j,-4)    = cell_u_tagB(i,j,1)
       cell_u_tagB(i,j,-5)    = cell_u_tagB(i,j,1)

       cell_u_tagB(i,j,k1)    = cell_u_tagB(i,j,kmax)
       cell_u_tagB(i,j,k1+1)  = cell_u_tagB(i,j,kmax)
       cell_u_tagB(i,j,k1+2)  = cell_u_tagB(i,j,kmax)
       cell_u_tagB(i,j,k1+3)  = cell_u_tagB(i,j,kmax)
       cell_u_tagB(i,j,k1+4)  = cell_u_tagB(i,j,kmax)
       cell_u_tagB(i,j,k1+5)  = cell_u_tagB(i,j,kmax)

       cell_v_tagB(i,j,0)     = cell_v_tagB(i,j,1)
       cell_v_tagB(i,j,-1)    = cell_v_tagB(i,j,1)
       cell_v_tagB(i,j,-2)    = cell_v_tagB(i,j,1)
       cell_v_tagB(i,j,-3)    = cell_v_tagB(i,j,1)
       cell_v_tagB(i,j,-4)    = cell_v_tagB(i,j,1)
       cell_v_tagB(i,j,-5)    = cell_v_tagB(i,j,1)

       cell_v_tagB(i,j,k1)    = cell_v_tagB(i,j,kmax)
       cell_v_tagB(i,j,k1+1)  = cell_v_tagB(i,j,kmax)
       cell_v_tagB(i,j,k1+2)  = cell_v_tagB(i,j,kmax)
       cell_v_tagB(i,j,k1+3)  = cell_v_tagB(i,j,kmax)
       cell_v_tagB(i,j,k1+4)  = cell_v_tagB(i,j,kmax)
       cell_v_tagB(i,j,k1+5)  = cell_v_tagB(i,j,kmax)

       cell_w_tagB(i,j,0)     = cell_w_tagB(i,j,1)
       cell_w_tagB(i,j,-1)    = cell_w_tagB(i,j,1)
       cell_w_tagB(i,j,-2)    = cell_w_tagB(i,j,1)
       cell_w_tagB(i,j,-3)    = cell_w_tagB(i,j,1)
       cell_w_tagB(i,j,-4)    = cell_w_tagB(i,j,1)
       cell_w_tagB(i,j,-5)    = cell_w_tagB(i,j,1)

       cell_w_tagB(i,j,k1)    = cell_w_tagB(i,j,kmax)
       cell_w_tagB(i,j,k1+1)  = cell_w_tagB(i,j,kmax)
       cell_w_tagB(i,j,k1+2)  = cell_w_tagB(i,j,kmax)
       cell_w_tagB(i,j,k1+3)  = cell_w_tagB(i,j,kmax)
       cell_w_tagB(i,j,k1+4)  = cell_w_tagB(i,j,kmax)
       cell_w_tagB(i,j,k1+5)  = cell_w_tagB(i,j,kmax)

       nx_surfB(i,j,0)     = nx_surfB(i,j,1)
       nx_surfB(i,j,-1)    = nx_surfB(i,j,1)
       nx_surfB(i,j,-2)    = nx_surfB(i,j,1)
       nx_surfB(i,j,-3)    = nx_surfB(i,j,1)
       nx_surfB(i,j,-4)    = nx_surfB(i,j,1)
       nx_surfB(i,j,-5)    = nx_surfB(i,j,1)

       nx_surfB(i,j,k1)    = nx_surfB(i,j,kmax)
       nx_surfB(i,j,k1+1)  = nx_surfB(i,j,kmax)
       nx_surfB(i,j,k1+2)  = nx_surfB(i,j,kmax)
       nx_surfB(i,j,k1+3)  = nx_surfB(i,j,kmax)
       nx_surfB(i,j,k1+4)  = nx_surfB(i,j,kmax)
       nx_surfB(i,j,k1+5)  = nx_surfB(i,j,kmax)

       ny_surfB(i,j,0)     = ny_surfB(i,j,1)
       ny_surfB(i,j,-1)    = ny_surfB(i,j,1)
       ny_surfB(i,j,-2)    = ny_surfB(i,j,1)
       ny_surfB(i,j,-3)    = ny_surfB(i,j,1)
       ny_surfB(i,j,-4)    = ny_surfB(i,j,1)
       ny_surfB(i,j,-5)    = ny_surfB(i,j,1)

       ny_surfB(i,j,k1)    = ny_surfB(i,j,kmax)
       ny_surfB(i,j,k1+1)  = ny_surfB(i,j,kmax)
       ny_surfB(i,j,k1+2)  = ny_surfB(i,j,kmax)
       ny_surfB(i,j,k1+3)  = ny_surfB(i,j,kmax)
       ny_surfB(i,j,k1+4)  = ny_surfB(i,j,kmax)
       ny_surfB(i,j,k1+5)  = ny_surfB(i,j,kmax)


       nz_surfB(i,j,0)     = nz_surfB(i,j,1)
       nz_surfB(i,j,-1)    = nz_surfB(i,j,1)
       nz_surfB(i,j,-2)    = nz_surfB(i,j,1)
       nz_surfB(i,j,-3)    = nz_surfB(i,j,1)
       nz_surfB(i,j,-4)    = nz_surfB(i,j,1)
       nz_surfB(i,j,-5)    = nz_surfB(i,j,1)

       nz_surfB(i,j,k1)    = nz_surfB(i,j,kmax)
       nz_surfB(i,j,k1+1)  = nz_surfB(i,j,kmax)
       nz_surfB(i,j,k1+2)  = nz_surfB(i,j,kmax)
       nz_surfB(i,j,k1+3)  = nz_surfB(i,j,kmax)
       nz_surfB(i,j,k1+4)  = nz_surfB(i,j,kmax)
       nz_surfB(i,j,k1+5)  = nz_surfB(i,j,kmax)



       nabs_surfB(i,j,0)     = nabs_surfB(i,j,1)
       nabs_surfB(i,j,-1)    = nabs_surfB(i,j,1)
       nabs_surfB(i,j,-2)    = nabs_surfB(i,j,1)
       nabs_surfB(i,j,-3)    = nabs_surfB(i,j,1)
       nabs_surfB(i,j,-4)    = nabs_surfB(i,j,1)
       nabs_surfB(i,j,-5)    = nabs_surfB(i,j,1)

       nabs_surfB(i,j,k1)    = nabs_surfB(i,j,kmax)
       nabs_surfB(i,j,k1+1)  = nabs_surfB(i,j,kmax)
       nabs_surfB(i,j,k1+2)  = nabs_surfB(i,j,kmax)
       nabs_surfB(i,j,k1+3)  = nabs_surfB(i,j,kmax)
       nabs_surfB(i,j,k1+4)  = nabs_surfB(i,j,kmax)
       nabs_surfB(i,j,k1+5)  = nz_surfB(i,j,kmax)



    enddo
 enddo
call updthalosBig(cell_phi_tagB,1)
call updthalosBig(cell_phi_tagB,2)
call updthalosBig(cell_u_tagB,1)
call updthalosBig(cell_u_tagB,1)
call updthalosBig(cell_v_tagB,2)
call updthalosBig(cell_v_tagB,2)
call updthalosBig(cell_w_tagB,1)
call updthalosBig(cell_w_tagB,2)
call updthalosBig(nx_surfB,1)
call updthalosBig(nx_surfB,2)
call updthalosBig(ny_surfB,1)
call updthalosBig(ny_surfB,2)
call updthalosBig(nz_surfB,1)
call updthalosBig(nz_surfB,1)
call updthalosBig(nabs_surfB,2)
call updthalosBig(nabs_surfB,2)

end subroutine boundloadIBM



subroutine updthalos(var,dir)
use mpi
use mod_param
use mod_common_mpi
implicit none
real, dimension(0:,0:,0:), intent(inout) :: var
integer, intent(in) :: dir
!
!  This subroutine updates the halos that store info
!  from the neighboring computational sub-domain
!
select case(dir)
case(1) ! x direction
call MPI_SENDRECV(var(1,0,0),1,xhalo,left,0,   &
                    var(i1,0,0),1,xhalo,right,0, &
                    comm_cart,status,error)
call MPI_SENDRECV(var(imax,0,0),1,xhalo,right,0, &
                    var(0,0,0),1,xhalo,left,0,     &
                    comm_cart,status,error)
!call MPI_IRECV(var(0,0,0),1,xhalo,left,1, &
!               comm_cart,requests(2),error)
!call MPI_IRECV(var(i1,0,0),1,xhalo,right,0, &
!               comm_cart,requests(1),error)
!call MPI_ISSEND(var(imax,0,0),1,xhalo,right,1, &
!               comm_cart,requests(4),error)
!call MPI_ISSEND(var(1,0,0),1,xhalo,left,0, &
!               comm_cart,requests(3),error)
!call MPI_WAITALL(4, requests, statuses, error)
case(2) ! y direction
  call MPI_SENDRECV(var(0,1,0),1,yhalo,front,0, &
                    var(0,j1,0),1,yhalo,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(0,jmax,0),1,yhalo,back,0, &
                    var(0,0,0),1,yhalo,front,0,   &
                    comm_cart,status,error)
!call MPI_IRECV(var(0,j1,0),1,yhalo,back,0, &
!               comm_cart,requests(1),error)
!call MPI_IRECV(var(0,0,0),1,yhalo,front,1, &
!               comm_cart,requests(2),error)
!call MPI_ISSEND(var(0,1,0),1,yhalo,front,0, &
!               comm_cart,requests(3),error)
!call MPI_ISSEND(var(0,jmax,0),1,yhalo,back,1, &
!               comm_cart,requests(4),error)
!call MPI_WAITALL(4, requests, statuses, error)
end select
!
return
end subroutine updthalos
!
!
subroutine updthalosBig(var,dir)
use mpi
use mod_param
use mod_common_mpi
implicit none
real, dimension(-5:i1+5,-5:j1+5,-5:k1+5), intent(inout) :: var
integer, intent(in) :: dir
!integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
!
!  This subroutine updates the halos that store info
!  from the neighboring computational sub-domain
!
select case(dir)
case(1) ! x direction

  call MPI_SENDRECV(var(1,-5,-5),1,xhalo3,left,0,   &
                    var(i1,-5,-5),1,xhalo3,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax,-5,-5),1,xhalo3,right,0, &
                    var(0,-5,-5),1,xhalo3,left,0,     &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(1+1,-5,-5),1,xhalo3,left,0,   &
                    var(i1+1,-5,-5),1,xhalo3,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax-1,-5,-5),1,xhalo3,right,0, &
                    var(0-1,-5,-5),1,xhalo3,left,0,     &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(1+2,-5,-5),1,xhalo3,left,0,   &
                    var(i1+2,-5,-5),1,xhalo3,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax-2,-5,-5),1,xhalo3,right,0, &
                    var(0-2,-5,-5),1,xhalo3,left,0,     &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(1+3,-5,-5),1,xhalo3,left,0,   &
                    var(i1+3,-5,-5),1,xhalo3,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax-3,-5,-5),1,xhalo3,right,0, &
                    var(0-3,-5,-5),1,xhalo3,left,0,     &
                    comm_cart,status,error)


  call MPI_SENDRECV(var(1+4,-5,-5),1,xhalo3,left,0,   &
                    var(i1+4,-5,-5),1,xhalo3,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax-4,-5,-5),1,xhalo3,right,0, &
                    var(0-4,-5,-5),1,xhalo3,left,0,     &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(1+5,-5,-5),1,xhalo3,left,0,   &
                    var(i1+5,-5,-5),1,xhalo3,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax-5,-5,-5),1,xhalo3,right,0, &
                    var(0-5,-5,-5),1,xhalo3,left,0,     &
                    comm_cart,status,error)


case(2) ! y direction

  call MPI_SENDRECV(var(-5,1,-5),1,yhalo3,front,0, &
                    var(-5,j1,-5),1,yhalo3,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-5,jmax,-5),1,yhalo3,back,0, &
                    var(-5,0,-5),1,yhalo3,front,0,   &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(-5,1+1,-5),1,yhalo3,front,0, &
                    var(-5,j1+1,-5),1,yhalo3,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-5,jmax-1,-5),1,yhalo3,back,0, &
                    var(-5,0-1,-5),1,yhalo3,front,0,   &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(-5,1+2,-5),1,yhalo3,front,0, &
                    var(-5,j1+2,-5),1,yhalo3,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-5,jmax-2,-5),1,yhalo3,back,0, &
                    var(-5,0-2,-5),1,yhalo3,front,0,   &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(-5,1+3,-5),1,yhalo3,front,0, &
                    var(-5,j1+3,-5),1,yhalo3,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-5,jmax-3,-5),1,yhalo3,back,0, &
                    var(-5,0-3,-5),1,yhalo3,front,0,   &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(-5,1+4,-5),1,yhalo3,front,0, &
                    var(-5,j1+4,-5),1,yhalo3,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-5,jmax-4,-5),1,yhalo3,back,0, &
                    var(-5,0-4,-5),1,yhalo3,front,0,   &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(-5,1+5,-5),1,yhalo3,front,0, &
                    var(-5,j1+5,-5),1,yhalo3,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-5,jmax-5,-5),1,yhalo3,back,0, &
                    var(-5,0-5,-5),1,yhalo3,front,0,   &
                    comm_cart,status,error)
end select
!
return
end subroutine updthalosBig

subroutine updthalosBigInt(var,dir)
use mpi
use mod_param
use mod_common_mpi
implicit none
integer, dimension(-5:i1+5,-5:j1+5,-5:k1+5), intent(inout) :: var
integer, intent(in) :: dir
!integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
!
!  This subroutine updates the halos that store info
!  from the neighboring computational sub-domain
!
select case(dir)
case(1) ! x direction

  call MPI_SENDRECV(var(1,-5,-5),1,xhalo2,left,0,   &
                    var(i1,-5,-5),1,xhalo2,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax,-5,-5),1,xhalo2,right,0, &
                    var(0,-5,-5),1,xhalo2,left,0,     &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(1+1,-5,-5),1,xhalo2,left,0,   &
                    var(i1+1,-5,-5),1,xhalo2,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax-1,-5,-5),1,xhalo2,right,0, &
                    var(0-1,-5,-5),1,xhalo2,left,0,     &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(1+2,-5,-5),1,xhalo2,left,0,   &
                    var(i1+2,-5,-5),1,xhalo2,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax-2,-5,-5),1,xhalo2,right,0, &
                    var(0-2,-5,-5),1,xhalo2,left,0,     &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(1+3,-5,-5),1,xhalo2,left,0,   &
                    var(i1+3,-5,-5),1,xhalo2,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax-3,-5,-5),1,xhalo2,right,0, &
                    var(0-3,-5,-5),1,xhalo2,left,0,     &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(1+4,-5,-5),1,xhalo2,left,0,   &
                    var(i1+4,-5,-5),1,xhalo2,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax-4,-5,-5),1,xhalo2,right,0, &
                    var(0-4,-5,-5),1,xhalo2,left,0,     &
                    comm_cart,status,error)


  call MPI_SENDRECV(var(1+5,-5,-5),1,xhalo2,left,0,   &
                    var(i1+5,-5,-5),1,xhalo2,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imax-5,-5,-5),1,xhalo2,right,0, &
                    var(0-5,-5,-5),1,xhalo2,left,0,     &
                    comm_cart,status,error)

case(2) ! y direction

  call MPI_SENDRECV(var(-5,1,-5),1,yhalo2,front,0, &
                    var(-5,j1,-5),1,yhalo2,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-5,jmax,-5),1,yhalo2,back,0, &
                    var(-5,0,-5),1,yhalo2,front,0,   &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(-5,1+1,-5),1,yhalo2,front,0, &
                    var(-5,j1+1,-5),1,yhalo2,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-5,jmax-1,-5),1,yhalo2,back,0, &
                    var(-5,0-1,-5),1,yhalo2,front,0,   &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(-5,1+2,-5),1,yhalo2,front,0, &
                    var(-5,j1+2,-5),1,yhalo2,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-5,jmax-2,-5),1,yhalo2,back,0, &
                    var(-5,0-2,-5),1,yhalo2,front,0,   &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(-5,1+3,-5),1,yhalo2,front,0, &
                    var(-5,j1+3,-5),1,yhalo2,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-5,jmax-3,-5),1,yhalo2,back,0, &
                    var(-5,0-3,-5),1,yhalo2,front,0,   &
                    comm_cart,status,error)


  call MPI_SENDRECV(var(-5,1+4,-5),1,yhalo2,front,0, &
                    var(-5,j1+4,-5),1,yhalo2,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-5,jmax-4,-5),1,yhalo2,back,0, &
                    var(-5,0-4,-5),1,yhalo2,front,0,   &
                    comm_cart,status,error)



  call MPI_SENDRECV(var(-5,1+5,-5),1,yhalo2,front,0, &
                    var(-5,j1+5,-5),1,yhalo2,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(-5,jmax-5,-5),1,yhalo2,back,0, &
                    var(-5,0-5,-5),1,yhalo2,front,0,   &
                    comm_cart,status,error)


end select


return
end subroutine updthalosBigInt



end module mod_bound
