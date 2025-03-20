module mod_output
  use mod_common
  use mod_common_mpi
  use mod_common_IBM
  use decomp_2d
  use decomp_2d_io
  use mod_post
  implicit none
  private
  public post3d
contains
  !
  !
  subroutine post3d(istep)
    implicit none
    integer, intent(in) :: istep
    real, allocatable, dimension(:,:,:) :: var1,var2
    integer, parameter :: nprocs=dims(1)*dims(2),ksol=kmax/nprocs
    integer :: p
    allocate(var1(1:imax,1:jmax,1:kmax))
    !
    ! pressure
    !
    !$omp workshare
    var1(:,:,:) = PFM_phi(1:imax,1:jmax,1:kmax)
    !$omp end workshare
    call write3dscal(istep,imax,jmax,kmax,var1,'phi',3,3,3)
    !
    ! u velocity
    !
    !$omp workshare
#ifdef IBM    
    var1(:,:,:) = cell_phi_tag(1:imax,1:jmax,1:kmax)
    !$omp end workshare
    call write3dscal(istep,imax,jmax,kmax,var1,'wal',3,3,3)
#endif
    ! v velocity
    !
    !$omp workshare
    var1(:,:,:) = 0.5*(vnew(1:imax,1:jmax,1:kmax)+vnew(1:imax,1-1:jmax-1,1:kmax))
    !$omp end workshare
    call write3dscal(istep,imax,jmax,kmax,var1,'vey',3,3,3)
    !
    ! w velocity
    !
    !$omp workshare
    var1(:,:,:) = 0.5*(wnew(1:imax,1:jmax,1:kmax)+wnew(1:imax,1:jmax,1-1:kmax-1))
    !$omp end workshare
    call write3dscal(istep,imax,jmax,kmax,var1,'vez',3,3,3)
    !
    ! dissipation (rate-of-strain)
    !
    !call strain_rate(unew,vnew,wnew,var1)
    !call write3dscal(istep,imax,jmax,kmax,var1,'str',4,4,4)
    !
    ! enstrophy
    !
    !call enstrophy(unew,vnew,wnew,var1)
    !call write3dscal(istep,imax,jmax,kmax,var1,'ens',4,4,4)
    !
    deallocate(var1)
    !allocate(var1(0:i1,0:j1,0:k1),var2(0:i1,0:j1,0:k1))
    !call phase_indicator(var2,var2,var2,var1, &
    !     var2,var2,var2,1)
    !
    ! particles' phase indicator
    !
    !call write3dscal(istep,imax,jmax,kmax,var1(1:imax,1:jmax,1:kmax),'sol',4,4,4)
    !call write3dscal(1,imax,jmax,kmax,unew(1:imax,1:jmax,1:kmax),'usm',4,4,4)
    !call write3dscal(1,imax,jmax,kmax,vnew(1:imax,1:jmax,1:kmax),'vsm',4,4,4)
    !call write3dscal(1,imax,jmax,kmax,wnew(1:imax,1:jmax,1:kmax),'wsm',4,4,4)
    !deallocate(var1,var2)
    !!
    !! Q-criterion
    !!
    !allocate(var3(1:imax,1:jmax,1:kmax))
    !call q_criterion(var2,var1,var3)
    !call write3dscal(istep,imax,jmax,kmax,var3,'qcr')
    !!
    !! R (third invariant of the velocity-gradient tensor
    !!
    !allocate(var4(1:imax,1:jmax,1:kmax))
    !call compute_r(unew,vnew,wnew,var4)
    !!
    !! swirling strength (lambda_{ci})
    !!
    !call swirl(var3,var4,var1)
    !call write3dscal(istep,imax,jmax,kmax,var3,'swr')
    !!
    return
  end subroutine post3d
  !
  subroutine write3dscal(istep,n1,n2,n3,var,name,iskip,jskip,kskip)
    implicit none
    integer, intent(in) :: istep,n1,n2,n3,iskip,jskip,kskip
    real, intent(in), dimension(n1,n2,n3) :: var
    character(len=3), intent(in) :: name
    integer :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    character :: istepchar*7
    !
    write(istepchar,'(i7.7)') istep
    !call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//name//istepchar, &
    !     MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    !filesize = 0_MPI_OFFSET_KIND
    !call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
    !disp = 0_MPI_OFFSET_KIND
    !call decomp_2d_write_var(fh,disp,3,var)
    !call MPI_FILE_CLOSE(fh,error)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//name//istepchar//'_3d.bin', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_write_every(3,var,iskip,jskip,kskip,datadir//name//istepchar//'_3d.bin',.true.)
    call MPI_FILE_CLOSE(fh,error)
    !
    return
  end subroutine write3dscal
end module mod_output

