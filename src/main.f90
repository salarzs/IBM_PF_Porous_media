program interpt
use mod_param
use mod_common
use mod_common_mpi
use mod_initmpi
use mod_init
use mod_bound
use mod_chkdiv
use mod_chkdt
use mod_loadflds
use mod_rk
use mod_initsolver
use mod_fillps
use mod_solver
use mod_correc
use mod_output
use mod_vtk_write
use mod_post
use mod_mom
use mod_interface
use mod_initIBM
use mod_IBM
use mod_common_IBM
use mod_debug
use mod_CahnHilliardSI
use mod_Impl_PredVel
implicit none
integer :: begin,nstep,istep
real :: dtmax,dt_scale
real, dimension(0:i1,0:j1,0:k1) :: durkold,dvrkold,dwrkold,dPFMold
real, dimension(0:i1,0:j1,0:k1) :: dPFM_couple_x,dPFM_couple_y,dPFM_couple_z
real, target, dimension(0:i1,0:j1,0:k1) :: p
real, pointer, dimension(:,:,:) :: ppp,ustar,vstar,wstar
integer :: i,j,k,l
real :: norm,norm_all,v_bulk_all,w_bulk_all
real :: distance,xxx,yyy,zzz,DropletsCenty,DropletsCentz
real,dimension(0:i1,0:j1,0:k1):: dPFM_visual, Phi_ghost
real:: time_star
logical::output1,output2,output3,output4,outputall
logical:: cond1, cond2, cond3,cond4,cond5,cond6
real::time_start,time_end

call initmpi
if (myid.eq.0)  write(6,*) '***************************************************'
if (myid.eq.0)  write(6,*) '***************************************************'
if (myid.eq.0)  write(6,*) 'Phase Field Code + Immersed Boundary Method '
if (myid.eq.0)  write(6,*) 'Surface tension Coefficient   = ', PFM_sigma_ref
if (myid.eq.0)  write(6,*) 'Viscosity                     = ', visc
if (myid.eq.0)  write(6,*) 'Interface thickness           = ', PFM_l
if (myid.eq.0)  write(6,*) 'Wall friction coefficient     = ', PFM_mu_f
if (myid.eq.0)  write(6,*) 'Gravity (Z direction)         = ', g_z
if (myid.eq.0)  write(6,*) 'Mobility coefficient          = ', mobility
if (myid.eq.0)  write(6,*) '***************************************************'
if (myid.eq.0)  write(6,*) '***************************************************'
#ifdef EXP
dt_scale=1e-2
#endif
#ifdef SI
dt_scale = 0.2
#endif
begin = 0 
nstep = 8000000 
!nstep = 4000000000 

#ifdef IBM
  call initIBM(begin)
#endif


if(begin.eq.0) then
  time = 0.
  call init
  call CPU_TIME(time_end)
  if (myid.eq.0) print*, 'Initialization time = ', (time_end-time_start),'S'
else
  unew(:,:,:) = 0.
  vnew(:,:,:) = 0.
  wnew(:,:,:) = 0.
  pnew(:,:,:) = 0.
  PFM_phi(:,:,:) = 0.
  psi(:,:,:)     =0.
  call loadflds(0,begin)
  call boundloadd(unew,vnew,wnew,pnew,PFM_phi)
  if (myid .eq. 0) write(6,*) 'nr steps at beginning simulation = ',begin
endif


call initsolver
ppp => p(1:imax,1:jmax,1:kmax) 
ustar => dudt(1:imax,1:jmax,1:kmax)
vstar => dvdt(1:imax,1:jmax,1:kmax)
wstar => dwdt(1:imax,1:jmax,1:kmax)
call filtering_phi
call boundPFM(PFM_phi,dPFM_bound_old,dPFM_bound_old_top,unew,vnew,wnew)
call filtering_phi
call ChemicalPotential
call boundChem(chem_pot)
call bounduvw(unew,vnew,wnew,PFM_phi,1)
call boundp(pnew)
call chkdiv
call chkdt(dtmax)
dt=1e-5 !dt_scale*dtmax
if (myid .eq. 0) write(6,*) 'dtmax = ', dtmax, ' dt = ', dt
durkold(:,:,:) = 0.
dvrkold(:,:,:) = 0.
dwrkold(:,:,:) = 0.
dPFMold(:,:,:) = 0.
dPFM_couple_x(:,:,:) = 0.
dPFM_couple_y(:,:,:) = 0.
dPFM_couple_z(:,:,:) = 0.
dPFM_bound_old(:,:)=0.0
dPFM_bound_old_top(:,:)=0.0
#ifdef IBM
dPFM_boundIBM(:,:,:) = 0.0
dPFM_boundIBM_old(:,:,:) = 0.0
#endif
! main loop below
do istep = begin+1,nstep
  time = time + dt
  time_star = time * PFM_omega
  u_bulk =sum(unew(1:imax,1:jmax,1:kmax))
  v_bulk =sum(vnew(1:imax,1:jmax,1:kmax))
  w_bulk =sum(wnew(1:imax,1:jmax,1:kmax))
  call mpi_allreduce(MPI_IN_PLACE,u_bulk,1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,v_bulk,1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(MPI_IN_PLACE,w_bulk,1,mpi_real8,mpi_sum,comm_cart,error)
  u_bulk=u_bulk/(1.*itot*jtot*ktot)
  v_bulk=v_bulk/(1.*itot*jtot*ktot)
  w_bulk=w_bulk/(1.*itot*jtot*ktot)
  if (myid.eq.0) write(6,*) 'time = ',time,'istep = ',istep,'ubulk=   ', u_bulk, 'vbulk = ',v_bulk,'wbulk = ',w_bulk
  poo        =   po  
  po         =   pnew
  uoo        =   uo  
  uo         =   unew
  voo        =   vo  
  vo         =   vnew
  woo        =   wo  
  wo         =   wnew 
  PFM_phioo  =   PFM_phio
  PFM_phio   =   PFM_phi
  if (istep.eq.begin+1) then
    poo        =  pnew
    po         =  pnew
    uoo        =  unew
    uo         =  unew
    voo        =  vnew
    vo         =  vnew
    woo        =  wnew
    wo         =  wnew
    PFM_phioo  =  PFM_phi
    PFM_phio   =  PFM_phi
  endif


!***********************************************************************
#ifdef EXP
  call updatePfm(dPFMold)
  call filtering_phi
  call boundPFM(PFM_phi,dPFM_bound_old,dPFM_bound_old_top,unew,vnew,wnew)
  call filtering_phi
  call ChemicalPotential
  call boundChem(chem_pot)
  call ab2(durkold,dvrkold,dwrkold,dPFM_couple_x,dPFM_couple_y,dPFM_couple_z)
#endif
#ifdef SI
  call updatePfmSemiImplicit(dt)
  call filtering_phi
  call boundPFM(PFM_phi,dPFM_bound_old,dPFM_bound_old_top,unew,vnew,wnew)
  call filtering_phi
  if (abs(PFM_mu_f).gt.1e-12) then 
   call updatePfmDyn
   call filtering_phi
   call boundPFM(PFM_phi,dPFM_bound_old,dPFM_bound_old_top,unew,vnew,wnew)
   call filtering_phi
  endif
  call ChemicalPotential
  call boundChem(chem_pot)
  call ab2_Imp(durkold,dvrkold,dwrkold)
  call bounduvw(dudt,dvdt,dwdt,PFM_phi,0)
  call Impl_PredVel(ustar,vstar,wstar)

#endif
!****************************************************
  call bounduvw(dudt,dvdt,dwdt,PFM_phi,1)
  !-----------------------------
#ifdef TwoD
  ustar = 0.
  dudt  = 0.
  unew  = 0.
#endif
  do k=1,kmax
   do j=1,jmax
    do i=1,imax
      if (istep.eq.begin+1) then
        phat(i,j,k) = po(i,j,k)
      else
        phat(i,j,k) = 2.0*po(i,j,k)-poo(i,j,k)
      endif
    enddo
   enddo
  enddo
  call boundp(phat)
  call fillps(p)
  call initsolver
  call solver2d(ppp,0.)
  call boundp(p)
  call correc(p)
  !-----------------------------
#ifdef TwoD
  ustar = 0.
  dudt  = 0.
  unew  = 0.
#endif
  !------------------------
  call bounduvw(unew,vnew,wnew,PFM_phi,0)
  pnew(:,:,:) = p(:,:,:)
  call boundp(pnew)
  k=kmax/2 ! near bottom wall
  norm = sum(sum(pnew(:,:,k),2),1)
  call mpi_allreduce(norm,norm_all,1,mpi_real8,mpi_sum,comm_cart,error)
  norm = norm_all/(1.*itot*jtot)
  pnew(:,:,:) = pnew(:,:,:) - norm
  cond1 = time_star.gt.9.1445
  cond2 = time_star.lt.9.1455
  cond3 = time_star.gt.9.555
  cond4 = time_star.lt.9.605
  cond5 = cond1.and.cond2
  cond6 = cond3.and.cond4
  if ((istep.eq.(begin+1)).or.(mod(istep,2000).eq.0)) then
#ifdef TwoD
     call fillVisual
     call post2dme(istep,1,itot/2,unew,vnew,wnew,pnew,C_Po)
#else
     call fillVisual
     call post2dme(istep,1,itot/2,unew,vnew,wnew,pnew,C_Po)
     call post2dme(istep,2,jtot/2,unew,vnew,wnew,pnew,C_Po)
     call post3d(istep)
#endif
  endif 
  call write_time(istep,time)
  if ((mod(istep,1000).eq.0).or.(istep.eq.(begin+1))) then
    call  Wetting_radius(time)
    call volume_of_droplet(time)
    call PrintBulkVelocities(time,u_bulk,v_bulk,w_bulk)

    !    call PFM_Energy 
  endif
  if (mod(istep,250).eq.0) then
    call chkdiv
    call chkdt(dtmax)
    !dt=dt_scale*dtmax
    dt = 1e-5
    if (myid .eq. 0) write(6,*) 'dtmax = ',dtmax,'dt = ',dt
  endif
  
  if (mod(istep,2000).eq.0) then
    call loadflds(1,istep)
  endif
enddo

if(myid.eq.0) write(6,*) '*** Fim ***'
!
call decomp_2d_finalize
call MPI_FINALIZE(error)
!
stop
end program 
