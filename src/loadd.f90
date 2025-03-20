module mod_loadflds
use decomp_2d
use decomp_2d_io
use mod_param 
use mod_common_mpi
use mod_common
use mod_common_IBM
implicit none
private
public loadflds, loadIBM
contains
!
subroutine loadflds(in,nr)
implicit none
integer :: in,nr
integer :: fh
integer(kind=MPI_OFFSET_KIND) :: filesize,disp
character(len=9) :: istepchar
real, dimension(3) :: fldinfo
real, dimension(imax,jmax,kmax) :: temp
integer:: i,j
if (in.eq.0) then
  write(istepchar,'(i9.9)') nr
  call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fld'//istepchar, &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
  call MPI_FILE_CLOSE(fh,error)
  call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fld'//istepchar, &
       MPI_MODE_RDONLY, MPI_INFO_NULL,fh, error)
  disp = 0_MPI_OFFSET_KIND
  call decomp_2d_read_var(fh,disp,3,temp)
  unew(1:imax,1:jmax,1:kmax) = temp
  call decomp_2d_read_var(fh,disp,3,temp)
  vnew(1:imax,1:jmax,1:kmax) = temp
  call decomp_2d_read_var(fh,disp,3,temp)
  wnew(1:imax,1:jmax,1:kmax) = temp
  call decomp_2d_read_var(fh,disp,3,temp)
  pnew(1:imax,1:jmax,1:kmax) = temp
  call decomp_2d_read_var(fh,disp,3,temp)
  PFM_phi(1:imax,1:jmax,1:kmax) = temp
  call decomp_2d_read_var(fh,disp,3,temp)
  do i=1,imax
   do j=1,jmax
    PFM_phi(i,j,0)   = temp(i,j,1)
    PFM_phi(i,j,k1)  = temp(i,j,2)
   enddo
  enddo
  call decomp_2d_read_scalar(fh,disp,3,fldinfo)
  time = fldinfo(1)
  nr = int(fldinfo(2))
  dt = fldinfo(3)
  call MPI_FILE_CLOSE(fh,error)
endif
!
if (in.eq.1) then
  write(istepchar,'(i9.9)') nr
  fldinfo = (/time,1.*nr,dt/)
  call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fld'//istepchar, &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
  filesize = 0_MPI_OFFSET_KIND
  call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
  disp = 0_MPI_OFFSET_KIND
  temp = unew(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
  temp = vnew(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
  temp = wnew(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
  temp = pnew(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
  temp = PFM_phi(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
  temp = 0.
  do i=1,imax
   do j=1,jmax
    temp(i,j,1) = PFM_phi(i,j,0)
    temp(i,j,2) = PFM_phi(i,j,k1)
    enddo
  enddo
  call decomp_2d_write_var(fh,disp,3,temp)


  call decomp_2d_write_scalar(fh,disp,3,fldinfo)
  call MPI_FILE_CLOSE(fh,error)
endif
!
return
end subroutine loadflds




subroutine loadIBM(in)
implicit none
integer,intent(in):: in
integer :: fh
integer(kind=MPI_OFFSET_KIND) :: filesize,disp
real, dimension(imax,jmax,kmax) :: temp
if (in.eq.0) then
  call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'IBMRestart', &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
  call MPI_FILE_CLOSE(fh,error)
  call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'IBMRestart', &
       MPI_MODE_RDONLY, MPI_INFO_NULL,fh, error)
  disp = 0_MPI_OFFSET_KIND
!1---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  cell_phi_tag(1:imax,1:jmax,1:kmax) = temp
!2---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  cell_u_tag(1:imax,1:jmax,1:kmax) = temp
!3---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  cell_v_tag(1:imax,1:jmax,1:kmax) = temp
!4---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  cell_w_tag(1:imax,1:jmax,1:kmax) = temp
!5---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  level_set(1:imax,1:jmax,1:kmax) = int(temp)
!6---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  nx_surf(1:imax,1:jmax,1:kmax) = temp
!7---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  ny_surf(1:imax,1:jmax,1:kmax) = temp
!8---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  nz_surf(1:imax,1:jmax,1:kmax) = temp
!9---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  nabs_surf(1:imax,1:jmax,1:kmax) = temp
!10---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  deltan(1:imax,1:jmax,1:kmax) = temp
!11---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  i_IP1(1:imax,1:jmax,1:kmax) = int(temp)
!12---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  j_IP1(1:imax,1:jmax,1:kmax) = int(temp)
!13---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  k_IP1(1:imax,1:jmax,1:kmax) = int(temp)
!14---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  i_IP2(1:imax,1:jmax,1:kmax) = int(temp)
!15---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  j_IP2(1:imax,1:jmax,1:kmax) = int(temp)
!16---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  k_IP2(1:imax,1:jmax,1:kmax) = int(temp)
!17---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  WP1(1:imax,1:jmax,1:kmax,1) = temp
!18---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  WP1(1:imax,1:jmax,1:kmax,2) = temp
!19---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  WP1(1:imax,1:jmax,1:kmax,3) = temp
!20---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  WP1(1:imax,1:jmax,1:kmax,4) = temp
!21---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  WP1(1:imax,1:jmax,1:kmax,5) = temp
!22---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  WP1(1:imax,1:jmax,1:kmax,6) = temp
!23---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  WP1(1:imax,1:jmax,1:kmax,7) = temp
!24---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  WP2(1:imax,1:jmax,1:kmax,1) = temp
!25---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  WP2(1:imax,1:jmax,1:kmax,2) = temp
!26---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  WP2(1:imax,1:jmax,1:kmax,3) = temp
!27---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  WP2(1:imax,1:jmax,1:kmax,4) = temp
!28---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  WP2(1:imax,1:jmax,1:kmax,5) = temp
!29---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  WP2(1:imax,1:jmax,1:kmax,6) = temp
!30---------------------------------------------
  call decomp_2d_read_var(fh,disp,3,temp)
  WP2(1:imax,1:jmax,1:kmax,7) = temp
!---------------------------------------------

  call MPI_FILE_CLOSE(fh,error)
endif
!
if (in.eq.1) then
  call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'IBMRestart', &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
  filesize = 0_MPI_OFFSET_KIND
  call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
  disp = 0_MPI_OFFSET_KIND
!1----------------------------------------------
  temp = cell_phi_tag(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!2---------------------------------------------
  temp = cell_u_tag(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!3---------------------------------------------
  temp = cell_v_tag(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!4---------------------------------------------
  temp = cell_w_tag(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!5---------------------------------------------
  temp = 1.0000*level_set(1:imax,1:jmax,1:kmax) 
  call decomp_2d_write_var(fh,disp,3,temp)
!6---------------------------------------------
  temp = nx_surf(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!7---------------------------------------------
  temp = ny_surf(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!8---------------------------------------------
  temp = nz_surf(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!9---------------------------------------------
  temp = nabs_surf(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!10---------------------------------------------
  temp = deltan(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
!11---------------------------------------------
  temp = 1.0000*(i_IP1(1:imax,1:jmax,1:kmax))
  call decomp_2d_write_var(fh,disp,3,temp)
!12---------------------------------------------
  temp = 1.0000*(j_IP1(1:imax,1:jmax,1:kmax))
  call decomp_2d_write_var(fh,disp,3,temp)
!13---------------------------------------------
  temp = 1.0000*(k_IP1(1:imax,1:jmax,1:kmax))
  call decomp_2d_write_var(fh,disp,3,temp)
!14---------------------------------------------
  temp = 1.0000*(i_IP2(1:imax,1:jmax,1:kmax))
  call decomp_2d_write_var(fh,disp,3,temp)
!15---------------------------------------------
  temp = 1.0000*(j_IP2(1:imax,1:jmax,1:kmax))
  call decomp_2d_write_var(fh,disp,3,temp)
!16---------------------------------------------
  temp = 1.0000*(k_IP2(1:imax,1:jmax,1:kmax))
  call decomp_2d_write_var(fh,disp,3,temp)
!17---------------------------------------------
  temp = WP1(1:imax,1:jmax,1:kmax,1)
  call decomp_2d_write_var(fh,disp,3,temp)
!18---------------------------------------------
  temp = WP1(1:imax,1:jmax,1:kmax,2)
  call decomp_2d_write_var(fh,disp,3,temp)
!19--------------------------------------------
  temp = WP1(1:imax,1:jmax,1:kmax,3)
  call decomp_2d_write_var(fh,disp,3,temp)
!20---------------------------------------------
  temp = WP1(1:imax,1:jmax,1:kmax,4)
  call decomp_2d_write_var(fh,disp,3,temp)
!21---------------------------------------------
  temp = WP1(1:imax,1:jmax,1:kmax,5)
  call decomp_2d_write_var(fh,disp,3,temp)
!22---------------------------------------------
  temp = WP1(1:imax,1:jmax,1:kmax,6)
  call decomp_2d_write_var(fh,disp,3,temp)
!23---------------------------------------------
  temp = WP1(1:imax,1:jmax,1:kmax,7)
  call decomp_2d_write_var(fh,disp,3,temp)
!24---------------------------------------------
  temp = WP2(1:imax,1:jmax,1:kmax,1)
  call decomp_2d_write_var(fh,disp,3,temp)
!25---------------------------------------------
  temp = WP2(1:imax,1:jmax,1:kmax,2)
  call decomp_2d_write_var(fh,disp,3,temp)
!26---------------------------------------------
  temp = WP2(1:imax,1:jmax,1:kmax,3)
  call decomp_2d_write_var(fh,disp,3,temp)
!27---------------------------------------------
  temp = WP2(1:imax,1:jmax,1:kmax,4)
  call decomp_2d_write_var(fh,disp,3,temp)
!28---------------------------------------------
  temp = WP2(1:imax,1:jmax,1:kmax,5)
  call decomp_2d_write_var(fh,disp,3,temp)
!29---------------------------------------------
  temp = WP2(1:imax,1:jmax,1:kmax,6)
  call decomp_2d_write_var(fh,disp,3,temp)
!30---------------------------------------------
  temp = WP2(1:imax,1:jmax,1:kmax,7)
  call decomp_2d_write_var(fh,disp,3,temp)
!---------------------------------------------

  call MPI_FILE_CLOSE(fh,error)
endif
!
return
end subroutine loadIBM






























end module mod_loadflds
