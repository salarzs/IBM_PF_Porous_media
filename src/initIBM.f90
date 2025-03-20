module mod_initIBM
use mod_IBM
use mod_bound
use mod_common
use decomp_2d
use mod_param
use mod_common_mpi
use mod_common_IBM
use mod_debug
use mod_loadflds 
implicit none
private
public initIBM
contains
subroutine initIBM(begin)
implicit none
integer,intent(in):: begin
integer::i,j,k
C_Po(:,:,:,:) = 0.
call allocate_arrays(begin)


if (begin.eq.0) then
  cell_u_tag(-5:i1+5,-5:j1+5,-5:k1+5)    = 1.
  cell_v_tag(-5:i1+5,-5:j1+5,-5:k1+5)    = 1.
  cell_w_tag(-5:i1+5,-5:j1+5,-5:k1+5)    = 1.
  cell_phi_tag(-5:i1+5,-5:j1+5,-5:k1+5)  = 1.
  Level_set(-5:i1+5,-5:j1+5,-5:k1+5)     = 1.
  call IBM_mask
      call updthalosBig(cell_u_tag,1)
      call updthalosBig(cell_u_tag,2)
      call updthalosBig(cell_v_tag,1)
      call updthalosBig(cell_v_tag,2)
      call updthalosBig(cell_w_tag,1)
      call updthalosBig(cell_w_tag,2)
      call updthalosBig(cell_phi_tag,1)
      call updthalosBig(cell_phi_tag,2)
      call updthalosBigInt(Level_set,1)
      call updthalosBigInt(Level_set,2)
       cell_u_tag(:,:,0)      = cell_u_tag(:,:,1)
       cell_u_tag(:,:,-1)     = cell_u_tag(:,:,1)
       cell_u_tag(:,:,-2)     = cell_u_tag(:,:,1)
       cell_u_tag(:,:,-3)     = cell_u_tag(:,:,1)
       cell_u_tag(:,:,-4)     = cell_u_tag(:,:,1)
       cell_u_tag(:,:,-5)     = cell_u_tag(:,:,1)
       cell_u_tag(:,:,k1)     = cell_u_tag(:,:,kmax)
       cell_u_tag(:,:,k1+1)   = cell_u_tag(:,:,kmax)
       cell_u_tag(:,:,k1+2)   = cell_u_tag(:,:,kmax)
       cell_u_tag(:,:,k1+3)   = cell_u_tag(:,:,kmax)
       cell_u_tag(:,:,k1+4)   = cell_u_tag(:,:,kmax)
       cell_u_tag(:,:,k1+5)   = cell_u_tag(:,:,kmax)
  
       cell_v_tag(:,:,0)      = cell_v_tag(:,:,1)
       cell_v_tag(:,:,-1)     = cell_v_tag(:,:,1)
       cell_v_tag(:,:,-2)     = cell_v_tag(:,:,1)
       cell_v_tag(:,:,-3)     = cell_v_tag(:,:,1)
       cell_v_tag(:,:,-4)     = cell_v_tag(:,:,1)
       cell_v_tag(:,:,-5)     = cell_v_tag(:,:,1)
       cell_v_tag(:,:,k1)     = cell_v_tag(:,:,kmax)
       cell_v_tag(:,:,k1+1)   = cell_v_tag(:,:,kmax)
       cell_v_tag(:,:,k1+2)   = cell_v_tag(:,:,kmax)
       cell_v_tag(:,:,k1+3)   = cell_v_tag(:,:,kmax)
       cell_v_tag(:,:,k1+4)   = cell_v_tag(:,:,kmax)
       cell_v_tag(:,:,k1+5)   = cell_v_tag(:,:,kmax)
  
       cell_w_tag(:,:,0)      = cell_w_tag(:,:,1)
       cell_w_tag(:,:,-1)     = cell_w_tag(:,:,1)
       cell_w_tag(:,:,-2)     = cell_w_tag(:,:,1)
       cell_w_tag(:,:,-3)     = cell_w_tag(:,:,1)
       cell_w_tag(:,:,-4)     = cell_w_tag(:,:,1)
       cell_w_tag(:,:,-5)     = cell_w_tag(:,:,1)
  
       cell_w_tag(:,:,k1)     = cell_w_tag(:,:,kmax)
       cell_w_tag(:,:,k1+1)   = cell_w_tag(:,:,kmax)
       cell_w_tag(:,:,k1+2)   = cell_w_tag(:,:,kmax)
       cell_w_tag(:,:,k1+3)   = cell_w_tag(:,:,kmax)
       cell_w_tag(:,:,k1+4)   = cell_w_tag(:,:,kmax)
       cell_w_tag(:,:,k1+5)   = cell_w_tag(:,:,kmax)
  
       cell_phi_tag(:,:,0)    = cell_phi_tag(:,:,1)
       cell_phi_tag(:,:,-1)   = cell_phi_tag(:,:,1)
       cell_phi_tag(:,:,-2  ) = cell_phi_tag(:,:,1)
       cell_phi_tag(:,:,-3  ) = cell_phi_tag(:,:,1)
       cell_phi_tag(:,:,-4  ) = cell_phi_tag(:,:,1)
       cell_phi_tag(:,:,-5  ) = cell_phi_tag(:,:,1)
       cell_phi_tag(:,:,k1  ) = cell_phi_tag(:,:,kmax)
       cell_phi_tag(:,:,k1+1) = cell_phi_tag(:,:,kmax)
       cell_phi_tag(:,:,k1+2) = cell_phi_tag(:,:,kmax)
       cell_phi_tag(:,:,k1+3) = cell_phi_tag(:,:,kmax)
       cell_phi_tag(:,:,k1+4) = cell_phi_tag(:,:,kmax)
       cell_phi_tag(:,:,k1+5) = cell_phi_tag(:,:,kmax)
  !     call MPI_WAIT(REQUEST,STATUS,error)
  if (myid.eq.0)  print*, 'Volume fractions have been calculated! '
  !---------------------------------------------------------------------
        Level_set(:,:,0)      = Level_set(:,:,1)
        Level_set(:,:,-1)     = Level_set(:,:,1)
        Level_set(:,:,-2)     = Level_set(:,:,1)
        Level_set(:,:,-3)     = Level_set(:,:,1)
        Level_set(:,:,-4)     = Level_set(:,:,1)
        Level_set(:,:,-5)     = Level_set(:,:,1)
        Level_set(:,:,k1)     = Level_set(:,:,kmax)
        Level_set(:,:,k1+1)   = Level_set(:,:,kmax)
        Level_set(:,:,k1+2)   = Level_set(:,:,kmax)
        Level_set(:,:,k1+3)   = Level_set(:,:,kmax)
        Level_set(:,:,k1+4)   = Level_set(:,:,kmax)
       Level_set(:,:,k1+5)   = Level_set(:,:,kmax)
  
  if (myid.eq.0)  print*, 'Level Set function has  been calculated! '
  !---------------------------------------------------------------------
  nx_surf(-5:i1+5,-5:j1+5,-5:k1+5)      = 0.
  ny_surf(-5:i1+5,-5:j1+5,-5:k1+5)      = 0.
  nz_surf(-5:i1+5,-5:j1+5,-5:k1+5)      = 0.
  nabs_surf(-5:i1+5,-5:j1+5,-5:k1+5)    = 0.
  call normal_vectors
  !call modify_normal_vectors
  
      call updthalosBig(nx_surf,1)
      call updthalosBig(nx_surf,2)
  
      call updthalosBig(ny_surf,1)
      call updthalosBig(ny_surf,2)
  
      call updthalosBig(nz_surf,1)
      call updthalosBig(nz_surf,2)
  
      call updthalosBig(nabs_surf,1)
      call updthalosBig(nabs_surf,2)
  
      nx_surf(:,:,0)          = nx_surf(:,:,1)
      nx_surf(:,:,-1)         = nx_surf(:,:,1)
      nx_surf(:,:,-2)         = nx_surf(:,:,1)
      nx_surf(:,:,-3)         = nx_surf(:,:,1)
      nx_surf(:,:,-4)         = nx_surf(:,:,1)
      nx_surf(:,:,-5)         = nx_surf(:,:,1)
      nx_surf(:,:,k1)         = nx_surf(:,:,kmax)
      nx_surf(:,:,k1+1)       = nx_surf(:,:,kmax)
      nx_surf(:,:,k1+2)       = nx_surf(:,:,kmax)
      nx_surf(:,:,k1+3)       = nx_surf(:,:,kmax)
      nx_surf(:,:,k1+4)       = nx_surf(:,:,kmax)
      nx_surf(:,:,k1+5)       = nx_surf(:,:,kmax)
  
      ny_surf(:,:,0)          = ny_surf(:,:,1)
      ny_surf(:,:,-1)         = ny_surf(:,:,1)
      ny_surf(:,:,-2)         = ny_surf(:,:,1)
      ny_surf(:,:,-3)         = ny_surf(:,:,1)
      ny_surf(:,:,-4)         = ny_surf(:,:,1)
      ny_surf(:,:,-5)         = ny_surf(:,:,1)
      ny_surf(:,:,k1)         = ny_surf(:,:,kmax)
      ny_surf(:,:,k1+1)       = ny_surf(:,:,kmax)
      ny_surf(:,:,k1+2)       = ny_surf(:,:,kmax)
      ny_surf(:,:,k1+3)       = ny_surf(:,:,kmax)
      ny_surf(:,:,k1+4)       = ny_surf(:,:,kmax)
      ny_surf(:,:,k1+5)       = ny_surf(:,:,kmax)
  
      nz_surf(:,:,0)          = nz_surf(:,:,1)
      nz_surf(:,:,-1)         = nz_surf(:,:,1)
      nz_surf(:,:,-2)         = nz_surf(:,:,1)
      nz_surf(:,:,-3)         = nz_surf(:,:,1)
      nz_surf(:,:,-4)         = nz_surf(:,:,1)
      nz_surf(:,:,-5)         = nz_surf(:,:,1)
      nz_surf(:,:,k1)         = nz_surf(:,:,kmax)
      nz_surf(:,:,k1+1)       = nz_surf(:,:,kmax)
      nz_surf(:,:,k1+2)       = nz_surf(:,:,kmax)
      nz_surf(:,:,k1+3)       = nz_surf(:,:,kmax)
      nz_surf(:,:,k1+4)       = nz_surf(:,:,kmax)
      nz_surf(:,:,k1+5)       = nz_surf(:,:,kmax)
  
      nabs_surf(:,:,0)        = nabs_surf(:,:,1)
      nabs_surf(:,:,-1)       = nabs_surf(:,:,1)
      nabs_surf(:,:,-2)       = nabs_surf(:,:,1)
      nabs_surf(:,:,-3)       = nabs_surf(:,:,1)
      nabs_surf(:,:,-4)       = nabs_surf(:,:,1)
      nabs_surf(:,:,-5)       = nabs_surf(:,:,1)
      nabs_surf(:,:,k1)       = nabs_surf(:,:,kmax)
      nabs_surf(:,:,k1+1)     = nabs_surf(:,:,kmax)
      nabs_surf(:,:,k1+2)     = nabs_surf(:,:,kmax)
      nabs_surf(:,:,k1+3)     = nabs_surf(:,:,kmax)
      nabs_surf(:,:,k1+4)     = nabs_surf(:,:,kmax)
      nabs_surf(:,:,k1+5)     = nabs_surf(:,:,kmax)
  if (myid.eq.0)  print*, 'Normal vectors have been calculated! '
  !call deallocate_volume_fractions
  !*********************************************************************
  x_intersect(-5:i1+5,-5:j1+5,-5:k1+5) = 0.
  y_intersect(-5:i1+5,-5:j1+5,-5:k1+5) = 0.
  z_intersect(-5:i1+5,-5:j1+5,-5:k1+5) = 0.
  call intersect
      call updthalosBig(x_intersect,1)
      call updthalosBig(x_intersect,2)
  
      call updthalosBig(y_intersect,1)
      call updthalosBig(y_intersect,2)
  
      call updthalosBig(z_intersect,1)
      call updthalosBig(z_intersect,2)
  
  
      x_intersect(:,:,0)        = x_intersect(:,:,1)
      x_intersect(:,:,-1)       = x_intersect(:,:,1)
      x_intersect(:,:,-2)       = x_intersect(:,:,1)
      x_intersect(:,:,-3)       = x_intersect(:,:,1)
      x_intersect(:,:,-4)       = x_intersect(:,:,1)
      x_intersect(:,:,-5)       = x_intersect(:,:,1)
      x_intersect(:,:,k1)       = x_intersect(:,:,kmax)
      x_intersect(:,:,k1+1)     = x_intersect(:,:,kmax)
      x_intersect(:,:,k1+2)     = x_intersect(:,:,kmax)
      x_intersect(:,:,k1+3)     = x_intersect(:,:,kmax)
      x_intersect(:,:,k1+4)     = x_intersect(:,:,kmax)
      x_intersect(:,:,k1+5)     = x_intersect(:,:,kmax)
  
      y_intersect(:,:,0)        = y_intersect(:,:,1)
      y_intersect(:,:,-1)       = y_intersect(:,:,1)
      y_intersect(:,:,-2)       = y_intersect(:,:,1)
      y_intersect(:,:,-3)       = y_intersect(:,:,1)
      y_intersect(:,:,-4)       = y_intersect(:,:,1)
      y_intersect(:,:,-5)       = y_intersect(:,:,1)
      y_intersect(:,:,k1)       = y_intersect(:,:,kmax)
      y_intersect(:,:,k1+1)     = y_intersect(:,:,kmax)
      y_intersect(:,:,k1+2)     = y_intersect(:,:,kmax)
      y_intersect(:,:,k1+3)     = y_intersect(:,:,kmax)
      y_intersect(:,:,k1+4)     = y_intersect(:,:,kmax)
      y_intersect(:,:,k1+5)     = y_intersect(:,:,kmax)
  
      z_intersect(:,:,0)        = z_intersect(:,:,1)
      z_intersect(:,:,-1)       = z_intersect(:,:,1)
      z_intersect(:,:,-2)       = z_intersect(:,:,1)
      z_intersect(:,:,-3)       = z_intersect(:,:,1)
      z_intersect(:,:,-4)       = z_intersect(:,:,1)
      z_intersect(:,:,-5)       = z_intersect(:,:,1)
      z_intersect(:,:,k1)       = z_intersect(:,:,kmax)
      z_intersect(:,:,k1+1)     = z_intersect(:,:,kmax)
      z_intersect(:,:,k1+2)     = z_intersect(:,:,kmax)
      z_intersect(:,:,k1+3)     = z_intersect(:,:,kmax)
      z_intersect(:,:,k1+4)     = z_intersect(:,:,kmax)
      z_intersect(:,:,k1+5)     = z_intersect(:,:,kmax)
  
  if (myid.eq.0)  print*, 'Intersect points have  been calculated! '
  !*********************************************************************
  x_mirror(-5:i1+5,-5:j1+5,-5:k1+5)    =   -1000.
  y_mirror(-5:i1+5,-5:j1+5,-5:k1+5)    =   -1000.
  z_mirror(-5:i1+5,-5:j1+5,-5:k1+5)    =   -1000.
  x_IP1(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000.
  x_IP2(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000.
  y_IP1(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000.
  y_IP2(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000.
  z_IP1(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000.
  z_IP2(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000.
  deltan(-5:i1+5,-5:j1+5,-5:k1+5)      =   -1000.
      call mirrorpoints
      
      call updthalosBig(deltan,1)
      call updthalosBig(deltan,2)
  
      call updthalosBig(x_mirror,1)
      call updthalosBig(x_mirror,2)
      call updthalosBig(y_mirror,1)
      call updthalosBig(y_mirror,2)
      call updthalosBig(z_mirror,1)
      call updthalosBig(z_mirror,2)
  
      call updthalosBig(x_IP1,1)
      call updthalosBig(x_IP1,2)
      call updthalosBig(x_IP2,1)
      call updthalosBig(x_IP2,2)
      call updthalosBig(y_IP1,1)
      call updthalosBig(y_IP1,2)
      call updthalosBig(y_IP2,1)
      call updthalosBig(y_IP2,2)
      call updthalosBig(z_IP1,1)
      call updthalosBig(z_IP1,2)
      call updthalosBig(z_IP2,1)
      call updthalosBig(z_IP2,2)
  
  
      deltan(:,:,0)           = deltan(:,:,1)
      deltan(:,:,-1)          = deltan(:,:,1)
      deltan(:,:,-2)          = deltan(:,:,1)
      deltan(:,:,-3)          = deltan(:,:,1)
      deltan(:,:,-4)          = deltan(:,:,1)
      deltan(:,:,-5)          = deltan(:,:,1)
      deltan(:,:,k1)          = deltan(:,:,kmax)
      deltan(:,:,k1+1)        = deltan(:,:,kmax)
      deltan(:,:,k1+2)        = deltan(:,:,kmax)
      deltan(:,:,k1+3)        = deltan(:,:,kmax)
      deltan(:,:,k1+4)        = deltan(:,:,kmax)
      deltan(:,:,k1+5)        = deltan(:,:,kmax)
  
  
  
  
      x_mirror(:,:,0)        = x_mirror(:,:,1)
      x_mirror(:,:,-1)       = x_mirror(:,:,1)
      x_mirror(:,:,-2)       = x_mirror(:,:,1)
      x_mirror(:,:,-3)       = x_mirror(:,:,1)
      x_mirror(:,:,-4)       = x_mirror(:,:,1)
      x_mirror(:,:,-5)       = x_mirror(:,:,1)
      x_mirror(:,:,k1)       = x_mirror(:,:,kmax)
      x_mirror(:,:,k1+1)     = x_mirror(:,:,kmax)
      x_mirror(:,:,k1+2)     = x_mirror(:,:,kmax)
      x_mirror(:,:,k1+3)     = x_mirror(:,:,kmax)
      x_mirror(:,:,k1+4)     = x_mirror(:,:,kmax)
      x_mirror(:,:,k1+5)     = x_mirror(:,:,kmax)
  
      x_IP1(:,:,0)           = x_IP1(:,:,1) 
      x_IP1(:,:,-1)          = x_IP1(:,:,1)
      x_IP1(:,:,-2)          = x_IP1(:,:,1)
      x_IP1(:,:,-3)          = x_IP1(:,:,1)
      x_IP1(:,:,-4)          = x_IP1(:,:,1)
      x_IP1(:,:,-5)          = x_IP1(:,:,1)
      x_IP1(:,:,k1)          = x_IP1(:,:,kmax)
      x_IP1(:,:,k1+1)        = x_IP1(:,:,kmax)
      x_IP1(:,:,k1+2)        = x_IP1(:,:,kmax)
      x_IP1(:,:,k1+3)        = x_IP1(:,:,kmax)
      x_IP1(:,:,k1+4)        = x_IP1(:,:,kmax)
      x_IP1(:,:,k1+5)        = x_IP1(:,:,kmax)
  
      x_IP2(:,:,0)           = x_IP2(:,:,1)
      x_IP2(:,:,-1)          = x_IP2(:,:,1)
      x_IP2(:,:,-2)          = x_IP2(:,:,1)
      x_IP2(:,:,-3)          = x_IP2(:,:,1)
      x_IP2(:,:,-4)          = x_IP2(:,:,1)
      x_IP2(:,:,-5)          = x_IP2(:,:,1)
      x_IP2(:,:,k1)          = x_IP2(:,:,kmax)
      x_IP2(:,:,k1+1)        = x_IP2(:,:,kmax)
      x_IP2(:,:,k1+2)        = x_IP2(:,:,kmax)
      x_IP2(:,:,k1+3)        = x_IP2(:,:,kmax)
      x_IP2(:,:,k1+4)        = x_IP2(:,:,kmax)
      x_IP2(:,:,k1+5)        = x_IP2(:,:,kmax)
  
      y_mirror(:,:,0)        = y_mirror(:,:,1)
      y_mirror(:,:,-1)       = y_mirror(:,:,1)
      y_mirror(:,:,-2)       = y_mirror(:,:,1)
      y_mirror(:,:,-3)       = y_mirror(:,:,1)
      y_mirror(:,:,-4)       = y_mirror(:,:,1)
      y_mirror(:,:,-5)       = y_mirror(:,:,1)
      y_mirror(:,:,k1)       = y_mirror(:,:,kmax)
      y_mirror(:,:,k1+1)     = y_mirror(:,:,kmax)
      y_mirror(:,:,k1+2)     = y_mirror(:,:,kmax)
      y_mirror(:,:,k1+3)     = y_mirror(:,:,kmax)
      y_mirror(:,:,k1+4)     = y_mirror(:,:,kmax)
      y_mirror(:,:,k1+5)     = y_mirror(:,:,kmax)
  
  
  
      y_IP1(:,:,0)           = y_IP1(:,:,1)
      y_IP1(:,:,-1)          = y_IP1(:,:,1)
      y_IP1(:,:,-2)          = y_IP1(:,:,1)
      y_IP1(:,:,-3)          = y_IP1(:,:,1)
      y_IP1(:,:,-4)          = y_IP1(:,:,1)
      y_IP1(:,:,-5)          = y_IP1(:,:,1)
      y_IP1(:,:,k1)          = y_IP1(:,:,kmax)
      y_IP1(:,:,k1+1)        = y_IP1(:,:,kmax)
      y_IP1(:,:,k1+2)        = y_IP1(:,:,kmax)
      y_IP1(:,:,k1+3)        = y_IP1(:,:,kmax)
      y_IP1(:,:,k1+4)        = y_IP1(:,:,kmax)
      y_IP1(:,:,k1+5)        = y_IP1(:,:,kmax)
  
      y_IP2(:,:,0)           = y_IP2(:,:,1)
      y_IP2(:,:,-1)          = y_IP2(:,:,1)
      y_IP2(:,:,-2)          = y_IP2(:,:,1)
      y_IP2(:,:,-3)          = y_IP2(:,:,1)
      y_IP2(:,:,-4)          = y_IP2(:,:,1)
      y_IP2(:,:,-5)          = y_IP2(:,:,1)
      y_IP2(:,:,k1)          = y_IP2(:,:,kmax)
      y_IP2(:,:,k1+1)        = y_IP2(:,:,kmax)
      y_IP2(:,:,k1+2)        = y_IP2(:,:,kmax)
      y_IP2(:,:,k1+3)        = y_IP2(:,:,kmax)
      y_IP2(:,:,k1+4)        = y_IP2(:,:,kmax)
      y_IP2(:,:,k1+5)        = y_IP2(:,:,kmax)
  
  
      z_mirror(:,:,0)        = z_mirror(:,:,1)
      z_mirror(:,:,-1)       = z_mirror(:,:,1)
      z_mirror(:,:,-2)       = z_mirror(:,:,1)
      z_mirror(:,:,-3)       = z_mirror(:,:,1)
      z_mirror(:,:,-4)       = z_mirror(:,:,1)
      z_mirror(:,:,-5)       = z_mirror(:,:,1)
      z_mirror(:,:,k1)       = z_mirror(:,:,kmax)
      z_mirror(:,:,k1+1)     = z_mirror(:,:,kmax)
      z_mirror(:,:,k1+2)     = z_mirror(:,:,kmax)
      z_mirror(:,:,k1+3)     = z_mirror(:,:,kmax)
      z_mirror(:,:,k1+4)     = z_mirror(:,:,kmax)
      z_mirror(:,:,k1+5)     = z_mirror(:,:,kmax)
  
      z_IP1(:,:,0)           = z_IP1(:,:,1)
      z_IP1(:,:,-1)          = z_IP1(:,:,1)
      z_IP1(:,:,-2)          = z_IP1(:,:,1)
      z_IP1(:,:,-3)          = z_IP1(:,:,1)
      z_IP1(:,:,-4)          = z_IP1(:,:,1)
      z_IP1(:,:,-5)          = z_IP1(:,:,1)
      z_IP1(:,:,k1)          = z_IP1(:,:,kmax)
      z_IP1(:,:,k1+1)        = z_IP1(:,:,kmax)
      z_IP1(:,:,k1+2)        = z_IP1(:,:,kmax)
      z_IP1(:,:,k1+3)        = z_IP1(:,:,kmax)
      z_IP1(:,:,k1+4)        = z_IP1(:,:,kmax)
      z_IP1(:,:,k1+5)        = z_IP1(:,:,kmax)
  
  
      z_IP2(:,:,0)           = z_IP2(:,:,1)
      z_IP2(:,:,-1)          = z_IP2(:,:,1)
      z_IP2(:,:,-2)          = z_IP2(:,:,1)
      z_IP2(:,:,-3)          = z_IP2(:,:,1)
      z_IP2(:,:,-4)          = z_IP2(:,:,1)
      z_IP2(:,:,-5)          = z_IP2(:,:,1)
      z_IP2(:,:,k1)          = z_IP2(:,:,kmax)
      z_IP2(:,:,k1+1)        = z_IP2(:,:,kmax)
      z_IP2(:,:,k1+2)        = z_IP2(:,:,kmax)
      z_IP2(:,:,k1+3)        = z_IP2(:,:,kmax)
      z_IP2(:,:,k1+4)        = z_IP2(:,:,kmax)
      z_IP2(:,:,k1+5)        = z_IP2(:,:,kmax)
  
  
  
  
  
     i_mirror(-5:i1+5,-5:j1+5,-5:k1+5)    =   -1000
     j_mirror(-5:i1+5,-5:j1+5,-5:k1+5)    =   -1000
     k_mirror(-5:i1+5,-5:j1+5,-5:k1+5)    =   -1000
     i_IP1(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000
     j_IP1(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000
     k_IP1(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000
     i_IP2(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000
     j_IP2(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000
     k_IP2(-5:i1+5,-5:j1+5,-5:k1+5)       =   -1000
  
  
  
  call mirrorpoints_ijk
      i_mirror(:,:,0)        = i_mirror(:,:,1)
      i_mirror(:,:,-1)       = i_mirror(:,:,1)
      i_mirror(:,:,-2)       = i_mirror(:,:,1)
      i_mirror(:,:,-3)       = i_mirror(:,:,1)
      i_mirror(:,:,-4)       = i_mirror(:,:,1)
      i_mirror(:,:,-5)       = i_mirror(:,:,1)
      i_mirror(:,:,k1)       = i_mirror(:,:,kmax)
      i_mirror(:,:,k1+1)     = i_mirror(:,:,kmax)
      i_mirror(:,:,k1+2)     = i_mirror(:,:,kmax)
      i_mirror(:,:,k1+3)     = i_mirror(:,:,kmax)
      i_mirror(:,:,k1+4)     = i_mirror(:,:,kmax)
      i_mirror(:,:,k1+5)     = i_mirror(:,:,kmax)
  
      i_IP1(:,:,0)           = i_IP1(:,:,1)
      i_IP1(:,:,-1)          = i_IP1(:,:,1)
      i_IP1(:,:,-2)          = i_IP1(:,:,1)
      i_IP1(:,:,-3)          = i_IP1(:,:,1)
      i_IP1(:,:,-4)          = i_IP1(:,:,1)
      i_IP1(:,:,-5)          = i_IP1(:,:,1)
      i_IP1(:,:,k1)          = i_IP1(:,:,kmax)
      i_IP1(:,:,k1+1)        = i_IP1(:,:,kmax)
      i_IP1(:,:,k1+2)        = i_IP1(:,:,kmax) 
      i_IP1(:,:,k1+3)        = i_IP1(:,:,kmax)
      i_IP1(:,:,k1+4)        = i_IP1(:,:,kmax)
      i_IP1(:,:,k1+5)        = i_IP1(:,:,kmax)
  
      i_IP2(:,:,0)           = i_IP2(:,:,1)
      i_IP2(:,:,-1)          = i_IP2(:,:,1)
      i_IP2(:,:,-2)          = i_IP2(:,:,1)
      i_IP2(:,:,-3)          = i_IP2(:,:,1)
      i_IP2(:,:,-4)          = i_IP2(:,:,1)
      i_IP2(:,:,-5)          = i_IP2(:,:,1)
      i_IP2(:,:,k1)          = i_IP2(:,:,kmax)
      i_IP2(:,:,k1+1)        = i_IP2(:,:,kmax)
      i_IP2(:,:,k1+2)        = i_IP2(:,:,kmax) 
      i_IP2(:,:,k1+3)        = i_IP2(:,:,kmax)
      i_IP2(:,:,k1+4)        = i_IP2(:,:,kmax)
      i_IP2(:,:,k1+5)        = i_IP2(:,:,kmax)
  
      j_mirror(:,:,0)        = j_mirror(:,:,1)
      j_mirror(:,:,-1)       = j_mirror(:,:,1)
      j_mirror(:,:,-2)       = j_mirror(:,:,1)
      j_mirror(:,:,-3)       = j_mirror(:,:,1)
      j_mirror(:,:,-4)       = j_mirror(:,:,1)
      j_mirror(:,:,-5)       = j_mirror(:,:,1)
      j_mirror(:,:,k1)       = j_mirror(:,:,kmax)
      j_mirror(:,:,k1+1)     = j_mirror(:,:,kmax)
      j_mirror(:,:,k1+2)     = j_mirror(:,:,kmax)
      j_mirror(:,:,k1+3)     = j_mirror(:,:,kmax)
      j_mirror(:,:,k1+4)     = j_mirror(:,:,kmax)
      j_mirror(:,:,k1+5)     = j_mirror(:,:,kmax)
  
  
      j_IP1(:,:,0)           = j_IP1(:,:,1)
      j_IP1(:,:,-1)          = j_IP1(:,:,1)
      j_IP1(:,:,-2)          = j_IP1(:,:,1)
      j_IP1(:,:,-3)          = j_IP1(:,:,1)
      j_IP1(:,:,-4)          = j_IP1(:,:,1)
      j_IP1(:,:,-5)          = j_IP1(:,:,1)
      j_IP1(:,:,k1)          = j_IP1(:,:,kmax)
      j_IP1(:,:,k1+1)        = j_IP1(:,:,kmax)
      j_IP1(:,:,k1+2)        = j_IP1(:,:,kmax)    
      j_IP1(:,:,k1+3)        = j_IP1(:,:,kmax)
      j_IP1(:,:,k1+4)        = j_IP1(:,:,kmax)
      j_IP1(:,:,k1+5)        = j_IP1(:,:,kmax)
  
      j_IP2(:,:,0)           = j_IP2(:,:,1)
      j_IP2(:,:,-1)          = j_IP2(:,:,1)
      j_IP2(:,:,-2)          = j_IP2(:,:,1)
      j_IP2(:,:,-3)          = j_IP2(:,:,1)
      j_IP2(:,:,-4)          = j_IP2(:,:,1)
      j_IP2(:,:,-5)          = j_IP2(:,:,1)
      j_IP2(:,:,k1)          = j_IP2(:,:,kmax)
      j_IP2(:,:,k1+1)        = j_IP2(:,:,kmax)
      j_IP2(:,:,k1+2)        = j_IP2(:,:,kmax)   
      j_IP2(:,:,k1+3)        = j_IP2(:,:,kmax)
      j_IP2(:,:,k1+4)        = j_IP2(:,:,kmax)
      j_IP2(:,:,k1+5)        = j_IP2(:,:,kmax)
  
      k_mirror(:,:,0)        = k_mirror(:,:,1)
      k_mirror(:,:,-1)       = k_mirror(:,:,1)
      k_mirror(:,:,-2)       = k_mirror(:,:,1)
      k_mirror(:,:,-3)       = k_mirror(:,:,1)
      k_mirror(:,:,-4)       = k_mirror(:,:,1)
      k_mirror(:,:,-5)       = k_mirror(:,:,1)
      k_mirror(:,:,k1)       = k_mirror(:,:,kmax)
      k_mirror(:,:,k1+1)     = k_mirror(:,:,kmax)
      k_mirror(:,:,k1+2)     = k_mirror(:,:,kmax)
      k_mirror(:,:,k1+3)     = k_mirror(:,:,kmax)
      k_mirror(:,:,k1+4)     = k_mirror(:,:,kmax)
      k_mirror(:,:,k1+5)     = k_mirror(:,:,kmax)
  
      k_IP1(:,:,0)           = k_IP1(:,:,1)
      k_IP1(:,:,-1)          = k_IP1(:,:,1)
      k_IP1(:,:,-2)          = k_IP1(:,:,1)
      k_IP1(:,:,-3)          = k_IP1(:,:,1)
      k_IP1(:,:,-4)          = k_IP1(:,:,1)
      k_IP1(:,:,-5)          = k_IP1(:,:,1)
      k_IP1(:,:,k1)          = k_IP1(:,:,kmax)
      k_IP1(:,:,k1+1)        = k_IP1(:,:,kmax)
      k_IP1(:,:,k1+2)        = k_IP1(:,:,kmax)    
      k_IP1(:,:,k1+3)        = k_IP1(:,:,kmax)
      k_IP1(:,:,k1+4)        = k_IP1(:,:,kmax)
      k_IP1(:,:,k1+5)        = k_IP1(:,:,kmax)
  
      k_IP2(:,:,0)           = k_IP2(:,:,1)
      k_IP2(:,:,-1)          = k_IP2(:,:,1)
      k_IP2(:,:,-2)          = k_IP2(:,:,1)
      k_IP2(:,:,-3)          = k_IP2(:,:,1)
      k_IP2(:,:,-4)          = k_IP2(:,:,1)
      k_IP2(:,:,-5)          = k_IP2(:,:,1)
      k_IP2(:,:,k1)          = k_IP2(:,:,kmax)
      k_IP2(:,:,k1+1)        = k_IP2(:,:,kmax)
      k_IP2(:,:,k1+2)        = k_IP2(:,:,kmax)   
      k_IP2(:,:,k1+3)        = k_IP2(:,:,kmax)
      k_IP2(:,:,k1+4)        = k_IP2(:,:,kmax)
      k_IP2(:,:,k1+5)        = k_IP2(:,:,kmax)
  
  
  
  !call PrintGhostMirror
  call deallocate_intersects
  call deallocate_mirror_points
  
  if (myid.eq.0)  print*, 'Mirror points have  been calculated! '
  !------------------------------------------------------------
  WP1(-5:i1+5,-5:j1+5,-5:k1+5,:) = 0.
  WP2(-5:i1+5,-5:j1+5,-5:k1+5,:) = 0.
  call InterpolationWeights
  if (myid.eq.0)  print*, 'Interpolation Weights have been calculated! '
  call deallocate_interpolation_points
  
  call loadIBM(1)
else
  call loadIBM(0)
      call updthalosBig(cell_u_tag,1)
      call updthalosBig(cell_u_tag,2)
      call updthalosBig(cell_v_tag,1)
      call updthalosBig(cell_v_tag,2)
      call updthalosBig(cell_w_tag,1)
      call updthalosBig(cell_w_tag,2)
      call updthalosBig(cell_phi_tag,1)
      call updthalosBig(cell_phi_tag,2)
      call updthalosBigInt(Level_set,1)
      call updthalosBigInt(Level_set,2)
       cell_u_tag(:,:,0)      = cell_u_tag(:,:,1)
       cell_u_tag(:,:,-1)     = cell_u_tag(:,:,1)
       cell_u_tag(:,:,-2)     = cell_u_tag(:,:,1)
       cell_u_tag(:,:,-3)     = cell_u_tag(:,:,1)
       cell_u_tag(:,:,-4)     = cell_u_tag(:,:,1)
       cell_u_tag(:,:,-5)     = cell_u_tag(:,:,1)
       cell_u_tag(:,:,k1)     = cell_u_tag(:,:,kmax)
       cell_u_tag(:,:,k1+1)   = cell_u_tag(:,:,kmax)
       cell_u_tag(:,:,k1+2)   = cell_u_tag(:,:,kmax)
       cell_u_tag(:,:,k1+3)   = cell_u_tag(:,:,kmax)
       cell_u_tag(:,:,k1+4)   = cell_u_tag(:,:,kmax)
       cell_u_tag(:,:,k1+5)   = cell_u_tag(:,:,kmax)

       cell_v_tag(:,:,0)      = cell_v_tag(:,:,1)
       cell_v_tag(:,:,-1)     = cell_v_tag(:,:,1)
       cell_v_tag(:,:,-2)     = cell_v_tag(:,:,1)
       cell_v_tag(:,:,-3)     = cell_v_tag(:,:,1)
       cell_v_tag(:,:,-4)     = cell_v_tag(:,:,1)
       cell_v_tag(:,:,-5)     = cell_v_tag(:,:,1)
       cell_v_tag(:,:,k1)     = cell_v_tag(:,:,kmax)
       cell_v_tag(:,:,k1+1)   = cell_v_tag(:,:,kmax)
       cell_v_tag(:,:,k1+2)   = cell_v_tag(:,:,kmax)
       cell_v_tag(:,:,k1+3)   = cell_v_tag(:,:,kmax)
       cell_v_tag(:,:,k1+4)   = cell_v_tag(:,:,kmax)
       cell_v_tag(:,:,k1+5)   = cell_v_tag(:,:,kmax)

       cell_w_tag(:,:,0)      = cell_w_tag(:,:,1)
       cell_w_tag(:,:,-1)     = cell_w_tag(:,:,1)
       cell_w_tag(:,:,-2)     = cell_w_tag(:,:,1)
       cell_w_tag(:,:,-3)     = cell_w_tag(:,:,1)
       cell_w_tag(:,:,-4)     = cell_w_tag(:,:,1)
       cell_w_tag(:,:,-5)     = cell_w_tag(:,:,1)

       cell_w_tag(:,:,k1)     = cell_w_tag(:,:,kmax)
       cell_w_tag(:,:,k1+1)   = cell_w_tag(:,:,kmax)
       cell_w_tag(:,:,k1+2)   = cell_w_tag(:,:,kmax)
       cell_w_tag(:,:,k1+3)   = cell_w_tag(:,:,kmax)
       cell_w_tag(:,:,k1+4)   = cell_w_tag(:,:,kmax)
       cell_w_tag(:,:,k1+5)   = cell_w_tag(:,:,kmax)

       cell_phi_tag(:,:,0)    = cell_phi_tag(:,:,1)
       cell_phi_tag(:,:,-1)   = cell_phi_tag(:,:,1)
       cell_phi_tag(:,:,-2  ) = cell_phi_tag(:,:,1)
       cell_phi_tag(:,:,-3  ) = cell_phi_tag(:,:,1)
       cell_phi_tag(:,:,-4  ) = cell_phi_tag(:,:,1)
       cell_phi_tag(:,:,-5  ) = cell_phi_tag(:,:,1)
       cell_phi_tag(:,:,k1  ) = cell_phi_tag(:,:,kmax)
       cell_phi_tag(:,:,k1+1) = cell_phi_tag(:,:,kmax)
       cell_phi_tag(:,:,k1+2) = cell_phi_tag(:,:,kmax)
       cell_phi_tag(:,:,k1+3) = cell_phi_tag(:,:,kmax)
       cell_phi_tag(:,:,k1+4) = cell_phi_tag(:,:,kmax)
       cell_phi_tag(:,:,k1+5) = cell_phi_tag(:,:,kmax)

       Level_set(:,:,0)      = Level_set(:,:,1)
       Level_set(:,:,-1)     = Level_set(:,:,1)
       Level_set(:,:,-2)     = Level_set(:,:,1)
       Level_set(:,:,-3)     = Level_set(:,:,1)
       Level_set(:,:,-4)     = Level_set(:,:,1)
       Level_set(:,:,-5)     = Level_set(:,:,1)
       Level_set(:,:,k1)     = Level_set(:,:,kmax)
       Level_set(:,:,k1+1)   = Level_set(:,:,kmax)
       Level_set(:,:,k1+2)   = Level_set(:,:,kmax)
       Level_set(:,:,k1+3)   = Level_set(:,:,kmax)
       Level_set(:,:,k1+4)   = Level_set(:,:,kmax)
       Level_set(:,:,k1+5)   = Level_set(:,:,kmax)

       nx_surf(:,:,0)          = nx_surf(:,:,1)
       nx_surf(:,:,-1)         = nx_surf(:,:,1)
       nx_surf(:,:,-2)         = nx_surf(:,:,1)
       nx_surf(:,:,-3)         = nx_surf(:,:,1)
       nx_surf(:,:,-4)         = nx_surf(:,:,1)
       nx_surf(:,:,-5)         = nx_surf(:,:,1)
       nx_surf(:,:,k1)         = nx_surf(:,:,kmax)
       nx_surf(:,:,k1+1)       = nx_surf(:,:,kmax)
       nx_surf(:,:,k1+2)       = nx_surf(:,:,kmax)
       nx_surf(:,:,k1+3)       = nx_surf(:,:,kmax)
       nx_surf(:,:,k1+4)       = nx_surf(:,:,kmax)
       nx_surf(:,:,k1+5)       = nx_surf(:,:,kmax)
 
       ny_surf(:,:,0)          = ny_surf(:,:,1)
       ny_surf(:,:,-1)         = ny_surf(:,:,1)
       ny_surf(:,:,-2)         = ny_surf(:,:,1)
       ny_surf(:,:,-3)         = ny_surf(:,:,1)
       ny_surf(:,:,-4)         = ny_surf(:,:,1)
       ny_surf(:,:,-5)         = ny_surf(:,:,1)
       ny_surf(:,:,k1)         = ny_surf(:,:,kmax)
       ny_surf(:,:,k1+1)       = ny_surf(:,:,kmax)
       ny_surf(:,:,k1+2)       = ny_surf(:,:,kmax)
       ny_surf(:,:,k1+3)       = ny_surf(:,:,kmax)
       ny_surf(:,:,k1+4)       = ny_surf(:,:,kmax)
       ny_surf(:,:,k1+5)       = ny_surf(:,:,kmax)
 
       nz_surf(:,:,0)          = nz_surf(:,:,1)
       nz_surf(:,:,-1)         = nz_surf(:,:,1)
       nz_surf(:,:,-2)         = nz_surf(:,:,1)
       nz_surf(:,:,-3)         = nz_surf(:,:,1)
       nz_surf(:,:,-4)         = nz_surf(:,:,1)
       nz_surf(:,:,-5)         = nz_surf(:,:,1)
       nz_surf(:,:,k1)         = nz_surf(:,:,kmax)
       nz_surf(:,:,k1+1)       = nz_surf(:,:,kmax)
       nz_surf(:,:,k1+2)       = nz_surf(:,:,kmax)
       nz_surf(:,:,k1+3)       = nz_surf(:,:,kmax)
       nz_surf(:,:,k1+4)       = nz_surf(:,:,kmax)
       nz_surf(:,:,k1+5)       = nz_surf(:,:,kmax)
 
       nabs_surf(:,:,0)        = nabs_surf(:,:,1)
       nabs_surf(:,:,-1)       = nabs_surf(:,:,1)
       nabs_surf(:,:,-2)       = nabs_surf(:,:,1)
       nabs_surf(:,:,-3)       = nabs_surf(:,:,1)
       nabs_surf(:,:,-4)       = nabs_surf(:,:,1)
       nabs_surf(:,:,-5)       = nabs_surf(:,:,1)
       nabs_surf(:,:,k1)       = nabs_surf(:,:,kmax)
       nabs_surf(:,:,k1+1)     = nabs_surf(:,:,kmax)
       nabs_surf(:,:,k1+2)     = nabs_surf(:,:,kmax)
       nabs_surf(:,:,k1+3)     = nabs_surf(:,:,kmax)
       nabs_surf(:,:,k1+4)     = nabs_surf(:,:,kmax)
       nabs_surf(:,:,k1+5)     = nabs_surf(:,:,kmax)
 
      deltan(:,:,0)           = deltan(:,:,1)
      deltan(:,:,-1)          = deltan(:,:,1)
      deltan(:,:,-2)          = deltan(:,:,1)
      deltan(:,:,-3)          = deltan(:,:,1)
      deltan(:,:,-4)          = deltan(:,:,1)
      deltan(:,:,-5)          = deltan(:,:,1)
      deltan(:,:,k1)          = deltan(:,:,kmax)
      deltan(:,:,k1+1)        = deltan(:,:,kmax)
      deltan(:,:,k1+2)        = deltan(:,:,kmax)
      deltan(:,:,k1+3)        = deltan(:,:,kmax)
      deltan(:,:,k1+4)        = deltan(:,:,kmax)
      deltan(:,:,k1+5)        = deltan(:,:,kmax)


      i_IP1(:,:,0)           = i_IP1(:,:,1)
      i_IP1(:,:,-1)          = i_IP1(:,:,1)
      i_IP1(:,:,-2)          = i_IP1(:,:,1)
      i_IP1(:,:,-3)          = i_IP1(:,:,1)
      i_IP1(:,:,-4)          = i_IP1(:,:,1)
      i_IP1(:,:,-5)          = i_IP1(:,:,1)
      i_IP1(:,:,k1)          = i_IP1(:,:,kmax)
      i_IP1(:,:,k1+1)        = i_IP1(:,:,kmax)
      i_IP1(:,:,k1+2)        = i_IP1(:,:,kmax)
      i_IP1(:,:,k1+3)        = i_IP1(:,:,kmax)
      i_IP1(:,:,k1+4)        = i_IP1(:,:,kmax)
      i_IP1(:,:,k1+5)        = i_IP1(:,:,kmax)

      i_IP2(:,:,0)           = i_IP2(:,:,1)
      i_IP2(:,:,-1)          = i_IP2(:,:,1)
      i_IP2(:,:,-2)          = i_IP2(:,:,1)
      i_IP2(:,:,-3)          = i_IP2(:,:,1)
      i_IP2(:,:,-4)          = i_IP2(:,:,1)
      i_IP2(:,:,-5)          = i_IP2(:,:,1)
      i_IP2(:,:,k1)          = i_IP2(:,:,kmax)
      i_IP2(:,:,k1+1)        = i_IP2(:,:,kmax)
      i_IP2(:,:,k1+2)        = i_IP2(:,:,kmax)
      i_IP2(:,:,k1+3)        = i_IP2(:,:,kmax)
      i_IP2(:,:,k1+4)        = i_IP2(:,:,kmax)
      i_IP2(:,:,k1+5)        = i_IP2(:,:,kmax)

      j_IP1(:,:,0)           = j_IP1(:,:,1)
      j_IP1(:,:,-1)          = j_IP1(:,:,1)
      j_IP1(:,:,-2)          = j_IP1(:,:,1)
      j_IP1(:,:,-3)          = j_IP1(:,:,1)
      j_IP1(:,:,-4)          = j_IP1(:,:,1)
      j_IP1(:,:,-5)          = j_IP1(:,:,1)
      j_IP1(:,:,k1)          = j_IP1(:,:,kmax)
      j_IP1(:,:,k1+1)        = j_IP1(:,:,kmax)
      j_IP1(:,:,k1+2)        = j_IP1(:,:,kmax)
      j_IP1(:,:,k1+3)        = j_IP1(:,:,kmax)
      j_IP1(:,:,k1+4)        = j_IP1(:,:,kmax)
      j_IP1(:,:,k1+5)        = j_IP1(:,:,kmax)

      j_IP2(:,:,0)           = j_IP2(:,:,1)
      j_IP2(:,:,-1)          = j_IP2(:,:,1)
      j_IP2(:,:,-2)          = j_IP2(:,:,1)
      j_IP2(:,:,-3)          = j_IP2(:,:,1)
      j_IP2(:,:,-4)          = j_IP2(:,:,1)
      j_IP2(:,:,-5)          = j_IP2(:,:,1)
      j_IP2(:,:,k1)          = j_IP2(:,:,kmax)
      j_IP2(:,:,k1+1)        = j_IP2(:,:,kmax)
      j_IP2(:,:,k1+2)        = j_IP2(:,:,kmax)
      j_IP2(:,:,k1+3)        = j_IP2(:,:,kmax)
      j_IP2(:,:,k1+4)        = j_IP2(:,:,kmax)
      j_IP2(:,:,k1+5)        = j_IP2(:,:,kmax)


      k_IP1(:,:,0)           = k_IP1(:,:,1)
      k_IP1(:,:,-1)          = k_IP1(:,:,1)
      k_IP1(:,:,-2)          = k_IP1(:,:,1)
      k_IP1(:,:,-3)          = k_IP1(:,:,1)
      k_IP1(:,:,-4)          = k_IP1(:,:,1)
      k_IP1(:,:,-5)          = k_IP1(:,:,1)
      k_IP1(:,:,k1)          = k_IP1(:,:,kmax)
      k_IP1(:,:,k1+1)        = k_IP1(:,:,kmax)
      k_IP1(:,:,k1+2)        = k_IP1(:,:,kmax)
      k_IP1(:,:,k1+3)        = k_IP1(:,:,kmax)
      k_IP1(:,:,k1+4)        = k_IP1(:,:,kmax)
      k_IP1(:,:,k1+5)        = k_IP1(:,:,kmax)

      k_IP2(:,:,0)           = k_IP2(:,:,1)
      k_IP2(:,:,-1)          = k_IP2(:,:,1)
      k_IP2(:,:,-2)          = k_IP2(:,:,1)
      k_IP2(:,:,-3)          = k_IP2(:,:,1)
      k_IP2(:,:,-4)          = k_IP2(:,:,1)
      k_IP2(:,:,-5)          = k_IP2(:,:,1)
      k_IP2(:,:,k1)          = k_IP2(:,:,kmax)
      k_IP2(:,:,k1+1)        = k_IP2(:,:,kmax)
      k_IP2(:,:,k1+2)        = k_IP2(:,:,kmax)
      k_IP2(:,:,k1+3)        = k_IP2(:,:,kmax)
      k_IP2(:,:,k1+4)        = k_IP2(:,:,kmax)
      k_IP2(:,:,k1+5)        = k_IP2(:,:,kmax)



endif



return
end subroutine initIBM

subroutine allocate_arrays(begin)
implicit none

integer,intent(in):: begin
if (begin.eq.0) then
  allocate (z_intersect(-5:i1+5,-5:j1+5,-5:k1+5),y_intersect(-5:i1+5,-5:j1+5,-5:k1+5),&
               x_intersect(-5:i1+5,-5:j1+5,-5:k1+5))
  allocate(i_mirror(-5:i1+5,-5:j1+5,-5:k1+5),j_mirror(-5:i1+5,-5:j1+5,-5:k1+5), &
              k_mirror(-5:i1+5,-5:j1+5,-5:k1+5))
  allocate(x_mirror(-5:i1+5,-5:j1+5,-5:k1+5), y_mirror(-5:i1+5,-5:j1+5,-5:k1+5), &
              z_mirror(-5:i1+5,-5:j1+5,-5:k1+5))
  allocate(x_IP1(-5:i1+5,-5:j1+5,-5:k1+5), y_IP1(-5:i1+5,-5:j1+5,-5:k1+5), &
              z_IP1(-5:i1+5,-5:j1+5,-5:k1+5))
  allocate(x_IP2(-5:i1+5,-5:j1+5,-5:k1+5), y_IP2(-5:i1+5,-5:j1+5,-5:k1+5), &
              z_IP2(-5:i1+5,-5:j1+5,-5:k1+5))
endif
allocate(cell_u_tag(-5:i1+5,-5:j1+5,-5:k1+5),cell_v_tag(-5:i1+5,-5:j1+5,-5:k1+5), &
            cell_w_tag(-5:i1+5,-5:j1+5,-5:k1+5),cell_phi_tag(-5:i1+5,-5:j1+5,-5:k1+5))
allocate(Level_set(-5:i1+5,-5:j1+5,-5:k1+5))
allocate(nx_surf(-5:i1+5,-5:j1+5,-5:k1+5),ny_surf(-5:i1+5,-5:j1+5,-5:k1+5), &
           nz_surf(-5:i1+5,-5:j1+5,-5:k1+5),nabs_surf(-5:i1+5,-5:j1+5,-5:k1+5))
allocate(deltan(-5:i1+5,-5:j1+5,-5:k1+5))
allocate(i_IP1(-5:i1+5,-5:j1+5,-5:k1+5),j_IP1(-5:i1+5,-5:j1+5,-5:k1+5), &
            k_IP1(-5:i1+5,-5:j1+5,-5:k1+5))
allocate(i_IP2(-5:i1+5,-5:j1+5,-5:k1+5),j_IP2(-5:i1+5,-5:j1+5,-5:k1+5), &
            k_IP2(-5:i1+5,-5:j1+5,-5:k1+5))
allocate(WP1(-5:i1+5,-5:j1+5,-5:k1+5,7),WP2(-5:i1+5,-5:j1+5,-5:k1+5,7)) 

end subroutine allocate_arrays



!subroutine deallocate_volume_fractions

!deallocate(cell_u_tag,cell_v_tag, &
!            cell_w_tag,cell_phi_tag)
!deallocate(Level_set,Level_set_u,&
!            Level_set_v,Level_set_w)
!deallocate(nx_surf,ny_surf, &
!            nz_surf,nabs_surf)
!deallocate(nx_surf_u,ny_surf_u, &
!            nz_surf_u,nabs_surf_u)
!deallocate(nx_surf_v,ny_surf_v, &
!            nz_surf_v,nabs_surf_v)
!deallocate(nx_surf_w,ny_surf_w, &
!            nz_surf_w,nabs_surf_w)
!end subroutine deallocate_volume_fractions


subroutine deallocate_intersects
deallocate (z_intersect,y_intersect,&
             x_intersect)
end subroutine deallocate_intersects


subroutine deallocate_mirror_points 
deallocate(i_mirror,j_mirror, &
            k_mirror)
deallocate(x_mirror, y_mirror, &
            z_mirror)
end subroutine deallocate_mirror_points

subroutine deallocate_interpolation_points

deallocate(x_IP1,x_IP2,y_IP1, &
           y_IP2,z_IP1,z_IP2)
end subroutine deallocate_interpolation_points



end module mod_initIBM
