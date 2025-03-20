module mod_CahnHilliardSI
use mod_mom
use mod_common
use mod_common_mpi
use mod_bound
use mod_solver
use mod_initsolver
use mod_common_IBM
implicit none
private
public updatePfmSemiImplicit
contains
!
subroutine updatePfmSemiImplicit(dt)
implicit none
integer i,j,k,n
real, intent(in) :: dt
real ,dimension(0:i1,0:j1,0:k1) ::u_ast,v_ast,w_ast,temp_phi
real ,dimension(-5:i1+5,-5:j1+5,-5:k1+5) ::phi_ast,phi_hat
real:: ss,lambda,gamma1,gamma0
real:: term1 , term2, term3,temp1,temp2,temp3
real::temp_ip,temp_jp,temp_kp,temp_im,temp_jm,temp_km,temp
real::term_ip,term_jp,term_kp,term_im,term_jm,term_km,term
real:: helm_1,helm_2
real, dimension(1:imax,1:jmax,1:kmax) :: ppp
u_ast    = 2*uo-uoo
#ifdef TwoD
u_ast    = 0.
#endif
v_ast    = 2.*vo-voo
w_ast    = 2.*wo-woo
phi_ast  = 2.*PFM_phio-PFM_phioo
phi_hat  = 2.*PFM_phio-0.5*PFM_phioo
gamma0   = 1.5 
gamma1   = mobility
lambda   = (3./(2.*sqrt(2.)))*PFM_sigma_ref*PFM_l
ss       = 1.01*PFM_l*PFM_l*sqrt((4.*gamma0)/(lambda*gamma1*dt))
temp1    = -ss/(2.*PFM_l*PFM_l)
temp2    = 4.*gamma0*PFM_l**4.
temp3    = lambda*gamma1*dt*ss*ss
alpha_PF = temp1*(1.+sqrt(1.-temp2/temp3))
helm_1   = -(alpha_PF+ss/(PFM_l*PFM_l))
helm_2   =  alpha_PF
 do k=1,kmax
  do j=1,jmax
   do i=1,imax
    temp_ip =  0.5*(phi_ast(i+1,j,k)+phi_ast(i,j,k))*u_ast(i,j,k) 
    temp_jp =  0.5*(phi_ast(i,j+1,k)+phi_ast(i,j,k))*v_ast(i,j,k)
    temp_kp =  0.5*(phi_ast(i,j,k+1)+phi_ast(i,j,k))*w_ast(i,j,k)
    temp_im =  0.5*(phi_ast(i-1,j,k)+phi_ast(i,j,k))*u_ast(i-1,j,k)  
    temp_jm =  0.5*(phi_ast(i,j-1,k)+phi_ast(i,j,k))*v_ast(i,j-1,k)
    temp_km =  0.5*(phi_ast(i,j,k-1)+phi_ast(i,j,k))*w_ast(i,j,k-1)
 
 
    term1   =  -(temp_ip-temp_im)*dxi-(temp_jp-temp_jm)*dyi-(temp_kp-temp_km)*dzi+phi_hat(i,j,k)/dt
 
    term_ip = (1./(PFM_l*PFM_l))*phi_ast(i+1,j,k)*(phi_ast(i+1,j,k)*phi_ast(i+1,j,k)-1.)- &
              ss/(PFM_l*PFM_l)*phi_ast(i+1,j,k)
    term_im = (1./(PFM_l*PFM_l))*phi_ast(i-1,j,k)*(phi_ast(i-1,j,k)*phi_ast(i-1,j,k)-1.)- &
               ss/(PFM_l*PFM_l)*phi_ast(i-1,j,k)
 
    term_jp = (1./(PFM_l*PFM_l))*phi_ast(i,j+1,k)*(phi_ast(i,j+1,k)*phi_ast(i,j+1,k)-1.)- &
               ss/(PFM_l*PFM_l)*phi_ast(i,j+1,k)
 
 
    term_jm = (1./(PFM_l*PFM_l))*phi_ast(i,j-1,k)*(phi_ast(i,j-1,k)*phi_ast(i,j-1,k)-1.)- &
              ss/(PFM_l*PFM_l)*phi_ast(i,j-1,k)
 
 
    term_kp = (1./(PFM_l*PFM_l))*phi_ast(i,j,k+1)*(phi_ast(i,j,k+1)*phi_ast(i,j,k+1)-1.)- &
              ss/(PFM_l*PFM_l)*phi_ast(i,j,k+1)
 
 
    term_km = (1./(PFM_l*PFM_l))*phi_ast(i,j,k-1)*(phi_ast(i,j,k-1)*phi_ast(i,j,k-1)-1.)- &
              ss/(PFM_l*PFM_l)*phi_ast(i,j,k-1)
 
 
    term = (1./(PFM_l*PFM_l))*phi_ast(i,j,k)*(phi_ast(i,j,k)*phi_ast(i,j,k)-1.)- &
               ss/(PFM_l*PFM_l)*phi_ast(i,j,k)
 
    term2 =  (term_ip-2.*term+term_im)*dxi*dxi+&
             (term_jp-2.*term+term_jm)*dyi*dyi+&
             (term_kp-2.*term+term_km)*dzi*dzi
    
    ppp(i,j,k) = (1./(lambda*gamma1))*term1+term2
 
   enddo
  enddo
 enddo
 
  call initsolver
  call solver2d(ppp,helm_1)
  
  call initsolver
  call solver2d(ppp,helm_2)

#ifdef IBM
    do k=1,kmax
     do j=1,jmax
      do i=1,imax
        if (nabs_surf(i,j,k).lt.1e-12) then
            PFM_phi(i,j,k) = ppp(i,j,k)
        endif
        if (cell_phi_tag(i,j,k).lt.1e-12) then
            PFM_phi(i,j,k) = -1.
        endif

      enddo
     enddo
    enddo
#endif
#ifdef NIBM
   PFM_phi(1:imax,1:jmax,1:kmax) = ppp(1:imax,1:jmax,1:kmax) 
#endif

return
end subroutine updatePfmSemiImplicit

end module mod_CahnHilliardSI
