module mod_IBM
use mod_param
use mod_common
use mod_common_mpi
use mod_common_IBM
implicit none
private
public IBM_Mask,normal_vectors,intersect,mirrorpoints,  &
        Penalization_center,Penalization_face,interpolation_2D_velocity,interpolation_2D_dphi, &
        mirrorpoints_ijk,modify_normal_vectors , interpolation_mirror,InterpolationWeights
contains
!

subroutine IBM_Mask
implicit none
integer i,j,k,l,n,m,number_of_divisions
real::xxx,yyy,zzz,dxx,dyy,dzz 
real:: cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z
real:: counter
logical :: inside,ghost
! Wall Geometry


number_of_divisions = 50 

do k=1,int(solid_height_ratio*k1)
  if (myid.eq.0) print*, 'Calculating volume fractions at k = ', k
  do j=1,jmax
    do i=1,imax
      xxx =  (i+coords(1)*imax)*dx-0.5*dx
      yyy =  (j+coords(2)*jmax)*dy-0.5*dy
      zzz =  (k)*dz - 0.5*dz
      call GuessGhostCells(xxx,yyy,zzz,ghost,inside)
      if (inside) then 
         cell_phi_tag(i,j,k) = 0.0
         cell_u_tag(i,j,k)   = 0.0
         cell_v_tag(i,j,k)   = 0.0
         cell_w_tag(i,j,k)   = 0.0
         Level_set(i,j,k)    = 0.0
      endif
      !if (.not.ghost)  cycle
      ! Cell Center
      inside = .false.
      cell_start_x = (i-1+coords(1)*imax)*dx
      cell_end_x = (i+coords(1)*imax)*dx

      cell_start_y = (j-1+coords(2)*jmax)*dy
      cell_end_y   = (j+coords(2)*jmax)*dy

      cell_start_z = (k-1)*dz
      cell_end_z   = (k)*dz

      dxx = (cell_end_x-cell_start_x)/number_of_divisions     
      dyy = (cell_end_y-cell_start_y)/number_of_divisions
      dzz = (cell_end_z-cell_start_z)/number_of_divisions

      counter = 0

      do n= 1,number_of_divisions
          zzz = cell_start_z+(n-1)*dzz
        do m = 1,number_of_divisions
            yyy = cell_start_y + (m-1)*dyy

            do l = 1,number_of_divisions
              xxx = cell_start_x + (l-1)*dxx
              call Solid_Surface(xxx,yyy,zzz,inside)
              if (inside) counter = counter +1
            enddo
        enddo
      enddo
     cell_phi_tag(i,j,k) = 1.0- counter/(1.0*number_of_divisions**3.)
! u cells
      inside = .false.

      cell_start_x = (i+coords(1)*imax)*dx-0.5*dx
      cell_end_x = (i+coords(1)*imax)*dx+0.5*dx

      cell_start_y = (j-1+coords(2)*jmax)*dy
      cell_end_y   = (j+coords(2)*jmax)*dy

      cell_start_z = (k-1)*dz
      cell_end_z   = (k)*dz
      dxx = (cell_end_x-cell_start_x)/number_of_divisions
      dyy = (cell_end_y-cell_start_y)/number_of_divisions
      dzz = (cell_end_z-cell_start_z)/number_of_divisions
      counter = 0 
      do n= 1,number_of_divisions
          xxx = 0 
          zzz = cell_start_z+(n-1)*dzz
        do m = 1,number_of_divisions
            yyy = cell_start_y + (m-1)*dyy
            do l = 1,number_of_divisions
              xxx = cell_start_x + (l-1)*dxx
              call Solid_Surface(xxx,yyy,zzz,inside)
              if (inside) counter = counter +1
            enddo
        enddo
      enddo

    cell_u_tag(i,j,k) = 1.0-counter/(number_of_divisions**3.)

! v cells
      inside = .false.

      cell_start_x = (i-1+coords(1)*imax)*dx
      cell_end_x = (i+coords(1)*imax)*dx

      cell_start_y = (j+coords(2)*jmax)*dy-0.5*dy
      cell_end_y   = (j+coords(2)*jmax)*dy+0.5*dy

      cell_start_z = (k-1)*dz
      cell_end_z   = (k)*dz
      dxx = (cell_end_x-cell_start_x)/number_of_divisions
      dyy = (cell_end_y-cell_start_y)/number_of_divisions
      dzz = (cell_end_z-cell_start_z)/number_of_divisions
      counter = 0 
      do n= 1,number_of_divisions
          xxx = 0 
          zzz = cell_start_z+(n-1)*dzz
        do m = 1,number_of_divisions
            yyy = cell_start_y + (m-1)*dyy
            do l = 1,number_of_divisions
              xxx = cell_start_x + (l-1)*dxx
              call Solid_Surface(xxx,yyy,zzz,inside)
              if (inside) counter = counter +1
            enddo
        enddo
      enddo

    cell_v_tag(i,j,k) = 1.0-counter/(number_of_divisions**3.)


! w cells
      inside = .false.

      cell_start_x = (i-1+coords(1)*imax)*dx
      cell_end_x = (i+coords(1)*imax)*dx

      cell_start_y = (j-1+coords(2)*jmax)*dy
      cell_end_y   = (j+coords(2)*jmax)*dy

      cell_start_z = (k)*dz-0.5*dz
      cell_end_z   = (k)*dz+0.5*dz
      dxx = (cell_end_x-cell_start_x)/number_of_divisions
      dyy = (cell_end_y-cell_start_y)/number_of_divisions
      dzz = (cell_end_z-cell_start_z)/number_of_divisions
      counter = 0 
      do n= 1,number_of_divisions
          xxx = 0 
          zzz = cell_start_z+(n-1)*dzz
        do m = 1,number_of_divisions
            yyy = cell_start_y + (m-1)*dyy
            do l = 1,number_of_divisions
              xxx = cell_start_x + (l-1)*dxx
              call Solid_Surface(xxx,yyy,zzz,inside)
            if (inside) counter = counter +1
           enddo
        enddo
      enddo

    cell_w_tag(i,j,k) = 1.0- counter/(number_of_divisions**3.)
    enddo
  enddo
enddo

end subroutine IBM_Mask




Subroutine GuessGhostCells(xxx,yyy,zzz,ghost,inside)
implicit none
real,intent(in)  :: xxx,yyy,zzz
logical,intent(out) :: ghost,inside
real:: RR, x_center, y_center, z_center,treshhold= 1.5*sqrt(dx*dx+dy*dy*dz*dz)
logical:: cond1,cond2,cond3,cond4,cond5
real:: LL, A, B
real::eps = 1.e-12
real:: x1,x2,x3,y1,y2,y3,y4,z1,z2,z3,z4
real::z_1,z_2,z_3,z_4,abscissa
real::y_begin, y_end,x_begin,x_end

select case(surface_type)


case('si3d')
     ghost=.false.
     inside= .false.
     y_begin = 0.05*ly
     y_end   = 0.95*ly
     x_begin = 0.05*lx
     x_end   = 0.95*lx
     inside= .false.
     A  = h_w*dz
     B  = h_b*dz
     x1 =  (lx/(2.*pi*n_sin))*asin((zzz/A)-(B/A)-sin(2.*pi*n_sin*yyy/ly))
     y1 =  (ly/(2.*pi*n_sin))*asin((zzz/A)-(B/A)-sin(2.*pi*n_sin*xxx/lx))
     z1 =  A*sin(n_sin*2*pi/ly*yyy)+A*sin(n_sin*2*pi/lx*xxx)+B
     cond1 = (abs(x1-xxx).lt.treshhold)
     cond2 = (abs(y1-yyy).lt.treshhold)
     cond3 = (abs(z1-zzz).lt.treshhold)
     if (cond1.or.cond2.or.cond3) ghost=.true.
     if ((yyy.lt.(y_begin)).or.(yyy.gt.(y_end)).or.(xxx.lt.(x_begin)).or.(xxx.gt.(x_end))) then
     LL = h_b*dz
     else
     LL = A*sin(n_sin*2*pi/ly*yyy)+A*sin(n_sin*2*pi/lx*xxx)+B
     endif
     if (zzz.le.LL) inside=.true.
case('sinu')
   inside= .false.
   LL = h_w*dz*sin(n_sin*2*pi/ly*yyy)+0.2*lz !h_b*dz
   if (zzz.le.LL) inside=.true.
   if (zzz.le.(0.5*lz)) ghost=.true. 

case('Sphe')
     ghost=.false.
     RR =  4*h_w*dz
     x_center = lx/2.
     y_center = ly/2.
     z_center = -0.75*RR
     x1 = x_center + sqrt(RR*RR-(zzz-z_center)*(zzz-z_center)-(yyy-y_center)*(yyy-y_center))
     x2 = x_center - sqrt(RR*RR-(zzz-z_center)*(zzz-z_center)-(yyy-y_center)*(yyy-y_center))

     y1 = y_center + sqrt(RR*RR-(zzz-z_center)*(zzz-z_center)-(xxx-x_center)*(xxx-x_center))
     y2 = y_center - sqrt(RR*RR-(zzz-z_center)*(zzz-z_center)-(xxx-x_center)*(xxx-x_center))

     z1 = z_center + sqrt(RR*RR-(xxx-x_center)*(xxx-x_center)-(yyy-y_center)*(yyy-y_center))
     z2 = z_center - sqrt(RR*RR-(xxx-x_center)*(xxx-x_center)-(yyy-y_center)*(yyy-y_center))
     cond1 = (abs(x1-xxx).lt.treshhold).or.(abs(x2-xxx).lt.treshhold)
     cond2 = (abs(y1-yyy).lt.treshhold).or.(abs(y2-yyy).lt.treshhold)
     cond3 = (abs(z1-zzz).lt.treshhold).or.(abs(z2-zzz).lt.treshhold)
     if (cond1.or.cond2.or.cond3) ghost=.true.

     inside=.false.
     LL = (zzz-z_center)*(zzz-z_center)+(yyy-y_center)*(yyy-y_center)+(xxx-x_center)*(xxx-x_center)
      if (LL.le.(RR*RR)) inside=.true.


end select



end Subroutine GuessGhostCells

Subroutine Solid_Surface(xxx,yyy,zzz,inside)

implicit none
real,intent(in)  :: xxx,yyy,zzz
logical,intent(out) :: inside
real:: LL, RR,x_begin,x_end, y_begin, y_end,x_center, y_center, z_center
real:: DD,L_periodic
integer::j_start,j_end,number_of_blocks ,l
real::z_1,z_2,z_3,z_4,b,abscissa,A
real:: x1,x2,y1,y2,y3,y4,z1,z2,z3,z4
logical:: cond1,cond2,cond3,cond4,cond5
real:: eps = 1e-15
real::y_a,z_a,y_b,z_b,y_c,z_c,y_d,z_d,y_e,z_e,y_f,z_f,y_g,z_g,y_h,z_h
real::yc1,zc1,yc2,zc2,yc3,zc3,yc4,zc4
real::z_5,z_6,z_7,z_8,z_9,z_10,z_11,z_12,z_13,z_14,z_15,z_16,z_17,z_18,z_19,z_20
logical :: smoothing= .false.
logical :: zone1,zone2,zone3,zone4
integer:: number_of_smoothed_grid = 8
real::y_cent_1,y_cent_2,y_cent_3,z_cent_1,z_cent_2,z_cent_3


inside= .false.

select case(surface_type)
case('circ')
   inside= .false.
   RR =  h_w*dz
   y_begin = ly/2 - RR
   y_end = y_begin+2*RR
   y_center = y_begin+ RR
   z_center = h_b*dz 
   
   LL = (zzz-z_center)*(zzz-z_center)+(yyy-y_center)*(yyy-y_center)
   if ((yyy.ge.y_begin).and.(yyy.le.y_end)) then
    if (LL.le.( (RR)*(RR))) inside=.true.
   endif 
    if (zzz.le.(h_b*dz)) inside=.true.

case('sinu')
   inside= .false.
   LL = h_w*dz*sin(n_sin*2*pi/ly*yyy)+0.2*lz!h_b*dz
   if (zzz.le.LL) inside=.true.


case('SiFl') !Flat plate with sinusoidal shape in the middle
   RR =  ly
   y_begin = ly/2 - 0.5*RR
   y_end = y_begin+RR
   
   inside= .false.
   if ((yyy.lt.(y_begin)).or.(yyy.gt.(y_end))) then
   LL = h_b*dz
   else
   LL = h_w*dz*sin(((2*n_sin-1)*pi/(RR))*(yyy- y_begin))+h_b*dz
   endif
   
   if (zzz.le.LL) inside=.true.
case('cosi')
   inside= .false.
   LL = h_w*dz*cos(n_sin*2*pi/ly*yyy)+h_b*dz
   !if (zzz.le.LL) inside=.true.
   if ((zzz-0.5*dz).le.LL) inside=.true.
case('CoFl') !Flat plate with sinusoidal shape in the middle
   RR =  ly  
   y_begin = 0.25*ly/n_sin
   y_end =  ly-y_begin
   inside= .false.
   if ((yyy.lt.(y_begin)).or.(yyy.gt.(y_end))) then
   LL = h_b*dz
   else
   LL = h_w*dz*cos(n_sin*2*pi/ly*(yyy+ y_begin))+h_b*dz 
   endif
   if (zzz.le.LL) inside=.true.

case('Flat')
   inside= .false.
   if (zzz.le.(h_w*dz)) inside=.true.


case('Poru')
inside = .false.
y_cent_1 = 0.25*ly
y_cent_2 = 0.50*ly
y_cent_3 = 0.75*ly

z_cent_1 = 0.25*lz
z_cent_2 = 0.50*lz
z_cent_3 = 0.75*lz

RR = 0.05*lz

cond1 = ((yyy-y_cent_1)**2+(zzz-z_cent_1)**2).LT.RR*RR
cond2 = ((yyy-y_cent_1)**2+(zzz-z_cent_3)**2).LT.RR*RR
cond3 = ((yyy-y_cent_2)**2+(zzz-z_cent_2)**2).LT.RR*RR
cond4 = ((yyy-y_cent_3)**2+(zzz-z_cent_1)**2).LT.RR*RR
cond5 = ((yyy-y_cent_3)**2+(zzz-z_cent_3)**2).LT.RR*RR




if (cond1.or.cond2.or.cond3.or.cond4.or.cond5) inside=.true.
case('RotB')
   inside= .false.
   A = 0.5*(ly-edge_length_of_box*(Cos_angle+Sin_angle))
   y1 = A + edge_length_of_box*Sin_angle
   z1 = A
   y2 = ly-A
   z2 = A + edge_length_of_box*Sin_angle
   y3 = A + edge_length_of_box*Cos_angle
   z3= lz - A
   y4 = A
   z4 = lz - (A+edge_length_of_box*Sin_angle)
   z_1 = ( (z2-z1)/(y2-y1+eps))*(yyy-y1)+z1
   z_2 = ( (z3-z2)/(y3-y2+eps))*(yyy-y2)+z2
   z_3 = ( (z4-z3)/(y4-y3+eps))*(yyy-y3)+z3
   z_4 = ( (z1-z4)/(y1-y4+eps))*(yyy-y4)+z4


   y_12_mid = 0.5*(y1+y2)
   z_12_mid = 0.5*(z1+z2)
   
   if (Rotation_angle.gt.1e-12) then 
    cond1 = zzz.gt.z_1
    cond2 = zzz.lt.z_3
    cond3 = zzz.gt.z_4
    cond4 = zzz.lt.z_2
   else
    cond1 = zzz.gt.z1
    cond2 = yyy.gt.y1
    cond3 = zzz.lt.z3
    cond4 = yyy.lt.y2
   endif
   cond5 = cond1.and.cond2.and.cond3.and.cond4
   if (.not.cond5) inside=.true.
   
   !-----------------------------------------------------------------
   ! Smoothing the corners
   if (smoothing) then
     RR = number_of_smoothed_grid*dz
     y_a = y1 + RR*cos_angle
     z_a = z1 + RR*sin_angle
     y_b = y2 - RR*cos_angle
     z_b = z2 - RR*sin_angle
     y_c = y2 - RR*sin_angle
     z_c = z2 + RR*cos_angle
     y_d = y3 + RR*sin_angle
     z_d = z3 - RR*cos_angle
     y_e = y3 - RR*cos_angle
     z_e = z3 - RR*sin_angle
     y_f = y4 + RR*cos_angle
     z_f = z4 + RR*sin_angle
     y_g = y4 + RR*sin_angle
     z_g = z4 - RR*cos_angle
     y_h = y1 - RR*sin_angle
     z_h = z1 + RR*cos_angle
     
     yc1 = y_h + RR*cos_angle
     zc1 = z_h + RR*sin_angle
     yc2 = y_c - RR*cos_angle
     zc2 = z_c - RR*sin_angle 
     yc3 = y_d - RR*cos_angle
     zc3 = z_d - RR*sin_angle
     yc4 = y_g + RR*cos_angle
     zc4 = z_g + RR*sin_angle
     
     z_5 = ( (zc1-z_h)/(yc1-y_h+eps))*(yyy-y_h)+z_h
     z_6 = ( (zc1-z_a)/(yc1-y_a+eps))*(yyy-y_a)+z_a
     z_7 = ( (zc2-z_b)/(yc2-y_b+eps))*(yyy-y_b)+z_b
     z_8 = ( (zc2-z_c)/(yc2-y_c+eps))*(yyy-y_c)+z_c
     z_9 = ( (zc3-z_d)/(yc3-y_d+eps))*(yyy-y_d)+z_d
     z_10 = ( (zc3-z_e)/(yc3-y_e+eps))*(yyy-y_e)+z_e
     z_11 = ( (zc4-z_f)/(yc4-y_f+eps))*(yyy-y_f)+z_f
     z_12 = ( (zc4-z_g)/(yc4-y_g+eps))*(yyy-y_g)+z_g
     
     z_13 = ( (z1-z_h)/(y1-y_h+eps))*(yyy-y_h)+z_h
     z_14 = ( (z1-z_a)/(y1-y_a+eps))*(yyy-y_a)+z_a
     z_15 = ( (z2-z_b)/(y2-y_b+eps))*(yyy-y_b)+z_b
     z_16 = ( (z2-z_c)/(y2-y_c+eps))*(yyy-y_c)+z_c
     z_17 = ( (z3-z_d)/(y3-y_d+eps))*(yyy-y_d)+z_d
     z_18 = ( (z3-z_e)/(y3-y_e+eps))*(yyy-y_e)+z_e
     z_19 = ( (z4-z_f)/(y4-y_f+eps))*(yyy-y_f)+z_f
     z_20 = ( (z4-z_g)/(y4-y_g+eps))*(yyy-y_g)+z_g
     if (Rotation_angle.gt.1e-12) then
     zone1 = (zzz.ge.z_14).and.(zzz.le.z_6).and.(zzz.ge.z_13).and.(zzz.le.z_5)
     zone2 = (zzz.ge.z_15).and.(zzz.le.z_16).and.(zzz.ge.z_7).and.(zzz.le.z_8)
     zone3 = (zzz.ge.z_9).and.(zzz.le.z_17).and.(zzz.ge.z_10).and.(zzz.le.z_18)
     zone4 = (zzz.ge.z_12).and.(zzz.le.z_11).and.(zzz.ge.z_20).and.(zzz.le.z_19)
     endif
     dd = 0.
     if (zone1) then
     inside = .false.
     dd =  (zzz-zc1)*(zzz-zc1)+(yyy-yc1)*(yyy-yc1)
     elseif (zone2) then
     dd =  (zzz-zc2)*(zzz-zc2)+(yyy-yc2)*(yyy-yc2)
     elseif (zone3) then
     dd =  (zzz-zc3)*(zzz-zc3)+(yyy-yc3)*(yyy-yc3)
     elseif (zone4) then
     dd =  (zzz-zc4)*(zzz-zc4)+(yyy-yc4)*(yyy-yc4)
     endif
     
     if (dd.gt.(RR*RR)) inside = .true.
   endif
case('Crev')
     inside= .false.
     abscissa = lz/2.
     b =edge_length_of_box *sqrt(2.)*0.4
     z1 =  (yyy - ly/2.) + lz/2.- b + abscissa
     z2 = -(yyy - ly/2.) + lz/2.- b + abscissa
     if (yyy.le.(0.5*ly)) then
          if ((zzz.le.z2)) then
              inside=.true.
          endif
     endif
     if (yyy.gt.(0.5*ly)) then
          if ((zzz.le.z1)) then
              inside=.true.
          endif
     endif

case('arcS')
     inside=.false.
     RR =  lz/3.!1.5
     y_begin = 0.!10 !ly/2 - RR
     y_end = ly !0.9 ! y_begin+2*RR
     y_center = ly/2!0.5 !y_begin+ RR
     z_center = 0.!-1.65+0.2+h_b*dz
     
     LL = (zzz-z_center)*(zzz-z_center)+(yyy-y_center)*(yyy-y_center)
     if ((yyy.gt.y_begin).and.(yyy.lt.y_end)) then
      if (LL.le.(RR*RR)) inside=.true.
     endif !else

case('jagg')
     inside=.false.
     LL = h_w*dz
     DD = h_b*dz
     L_periodic = (l_w+t_w)*dy
     number_of_blocks= int(ly/L_periodic)
     do l = 1, number_of_blocks+1
       j_start =  int(t_w/2) + (l-1)*(t_w+l_w)
       j_end   = l_w+int(t_w/2) + (l-1)*(t_w+l_w)
     
       y_begin = (j_start)*dy
       y_end   = (j_end)*dy
     
       if ((yyy.gt.y_begin).and.(yyy.lt.y_end)) then
          if (zzz.lt.LL) inside=.true.
       else
          if (zzz.lt.DD) inside=.true.
      endif
     enddo


case('Sphe')
     inside=.false.
     RR =  4*h_w*dz
     x_center = lx/2.
     y_center = ly/2.
     z_center = -0.75*RR

     LL = (zzz-z_center)*(zzz-z_center)+(yyy-y_center)*(yyy-y_center)+(xxx-x_center)*(xxx-x_center)
      if (LL.le.(RR*RR)) inside=.true.



case('si3d')
   inside= .false.
   y_begin = 0.1*ly
   y_end   = 0.9*ly
   x_begin = 0.1*lx
   x_end   = 0.9*lx
   inside= .false.
   if ((yyy.lt.(y_begin)).or.(yyy.gt.(y_end)).or.(xxx.lt.(x_begin)).or.(xxx.gt.(x_end))) then
   LL = h_b*dz
   else
     LL = h_w*dz*sin(n_sin*2*pi/ly*yyy)+h_w*dz*sin(n_sin*2*pi/lx*xxx)+h_b*dz
   endif
   if (zzz.le.LL) inside=.true.

end select

return

end subroutine Solid_Surface

!******************************************************************************************************************************************
Subroutine normal_vectors


implicit none
integer::i,j,k,l
real:: nx, ny, nz, n_abs,n_abs_p
real :: eps = 1e-12
real::m_av_x,m_av_y,m_av_z,normal_denum
real:: m1,m2,m3,m4,m5,m6,m7,m8
real:: cell_tag_wall = 0.5
integer:: ip,jp,kp
integer:: im,jm,km
integer:: iii,jjj,kkk
logical :: inside_1,inside_2,inside_3,inside_4,inside_5,inside_6,inside_7,inside_8,inside_9
logical :: inside_10,inside_11,inside_12,inside_13,inside_14,inside_15,inside_16,inside_17,inside_18
logical :: inside_19,inside_20,inside_21,inside_22,inside_23,inside_24,inside_25,inside_26,inside_27
logical ::inside_jpkm,inside_jmkp,inside_jmkm,inside_jk,ghost_cond, inside
real::xxx,yyy,zzz
integer, dimension(-5:i1+5,-5:j1+5,-5:k1+5) :: ghost_cell_tag 
integer:: Level_set_all
logical:: ghost
! Preparing normal vectors to the surfaces

! Identifying the ghost cells
ghost_cell_tag(:,:,:) = 0
do k=1,int(solid_height_ratio*k1)
if (myid.eq.0) print*, 'Calculating Normal Vectors at k = ', k
  do j=1,jmax
   do i=1,imax
            ghost = .false.
            xxx   =  (i+coords(1)*imax)*dx-0.5*dx
            yyy   =  (j+coords(2)*jmax)*dy-0.5*dy
            zzz   =  (k)*dz - 0.5*dz
            call GuessGhostCells(xxx,yyy,zzz,ghost,inside)
            if (.not.ghost)  cycle
            ip = i+1
            jp = j+1
            kp = k+1
            im = i-1
            jm = j-1
            km = k-1
            ghost_cond = .false.
            Level_set_all = Level_set(im,jm,kp) + Level_set(im,j,kp)+Level_set(im,jp,kp)  + & 
                       Level_set(i,jm,kp) + Level_set(i,j,kp)+Level_set(i,jp,kp) + & 
                       Level_set(ip,jm,kp) + Level_set(ip,j,kp)+Level_set(ip,jp,kp)+ &
                       Level_set(im,jm,k) + Level_set(im,j,k)+Level_set(im,jp,k)  + & 
                       Level_set(i,jm,k) + Level_set(i,j,k)+Level_set(i,jp,k) + & 
                       Level_set(ip,jm,k) + Level_set(ip,j,k)+Level_set(ip,jp,k)+ &
                       Level_set(im,jm,kp) + Level_set(im,j,kp)+Level_set(im,jp,kp)  + &   
                       Level_set(i,jm,km) + Level_set(i,j,km)+Level_set(i,jp,km) + &   
                       Level_set(ip,jm,km) + Level_set(ip,j,km)+Level_set(ip,jp,km)

            if ((Level_set_all.gt.0).and.(Level_set(i,j,k)).eq.0) ghost_cond =.true. 
            if(ghost_cond)  ghost_cell_tag(i,j,k) = 1
   enddo
 enddo
enddo

do k=1,kmax
 do j=1,jmax
  do i=1,imax
   if (ghost_cell_tag(i,j,k) .eq. 1) then
      nx = 0.
      ny = 0.
      nz = 0.
      n_abs= 0
      m_av_x = 0.
      m_av_y = 0.
      m_av_z = 0.
      normal_denum = 0.



      m1 =  (cell_phi_tag(i+1,j,k)+cell_phi_tag(i+1,j,k+1) &
            -cell_phi_tag(i,j,k)-cell_phi_tag(i,j,k+1)+ &
             cell_phi_tag(i+1,j+1,k)+cell_phi_tag(i+1,j+1,k+1) &
            -cell_phi_tag(i,j+1,k)-cell_phi_tag(i,j+1,k+1))*0.25*dxi


      m2 =  (cell_phi_tag(i,j,k)+cell_phi_tag(i,j,k+1) &
            -cell_phi_tag(i-1,j,k)-cell_phi_tag(i-1,j,k+1)+ &
             cell_phi_tag(i,j+1,k)+cell_phi_tag(i,j+1,k+1) &
            -cell_phi_tag(i-1,j+1,k)-cell_phi_tag(i-1,j+1,k+1))*0.25*dxi


      m3 =  (cell_phi_tag(i,j,k-1)+cell_phi_tag(i,j,k) &
            -cell_phi_tag(i-1,j,k-1)-cell_phi_tag(i-1,j,k)+ &
             cell_phi_tag(i,j+1,k-1)+cell_phi_tag(i,j+1,k) &
            -cell_phi_tag(i-1,j+1,k-1)-cell_phi_tag(i-1,j+1,k))*0.25*dxi

      m4 =  (cell_phi_tag(i+1,j,k-1)+cell_phi_tag(i+1,j,k) &
            -cell_phi_tag(i,j,k-1)-cell_phi_tag(i,j,k)+ &
             cell_phi_tag(i+1,j+1,k-1)+cell_phi_tag(i+1,j+1,k) &
            -cell_phi_tag(i,j+1,k-1)-cell_phi_tag(i,j+1,k))*0.25*dxi



      m5 =  (cell_phi_tag(i+1,j,k)+cell_phi_tag(i+1,j,k+1) &
            -cell_phi_tag(i,j,k)-cell_phi_tag(i,j,k+1)+ &
             cell_phi_tag(i+1,j-1,k)+cell_phi_tag(i+1,j-1,k+1) &
            -cell_phi_tag(i,j-1,k)-cell_phi_tag(i,j-1,k+1))*0.25*dxi


      m6 =  (cell_phi_tag(i,j,k)+cell_phi_tag(i,j,k+1) &
            -cell_phi_tag(i-1,j,k)-cell_phi_tag(i-1,j,k+1)+ &
             cell_phi_tag(i,j-1,k)+cell_phi_tag(i,j-1,k+1) &
            -cell_phi_tag(i-1,j-1,k)-cell_phi_tag(i-1,j-1,k+1))*0.25*dxi


      m7 =  (cell_phi_tag(i,j,k-1)+cell_phi_tag(i,j,k) &
            -cell_phi_tag(i-1,j,k-1)-cell_phi_tag(i-1,j,k)+ &
             cell_phi_tag(i,j-1,k-1)+cell_phi_tag(i,j-1,k) &
            -cell_phi_tag(i-1,j-1,k-1)-cell_phi_tag(i-1,j-1,k))*0.25*dxi

      m8 =  (cell_phi_tag(i+1,j,k-1)+cell_phi_tag(i+1,j,k) &
            -cell_phi_tag(i,j,k-1)-cell_phi_tag(i,j,k)+ &
             cell_phi_tag(i+1,j-1,k-1)+cell_phi_tag(i+1,j-1,k) &
            -cell_phi_tag(i,j-1,k-1)-cell_phi_tag(i,j-1,k))*0.25*dxi

      m_av_x= 0.125*(m1+m2+m3+m4+m5+m6+m7+m8)

      m1 = 0
      m2 = 0
      m3 = 0
      m4 = 0
      m5 = 0
      m6 = 0
      m7 = 0
      m8 = 0

      m1 =  (cell_phi_tag(i,j+1,k)+cell_phi_tag(i,j+1,k+1) &
            -cell_phi_tag(i,j,k)-cell_phi_tag(i,j,k+1)+ &
             cell_phi_tag(i+1,j+1,k)+cell_phi_tag(i+1,j+1,k+1) &
            -cell_phi_tag(i+1,j,k)-cell_phi_tag(i+1,j,k+1))*0.25*dyi


      m2 =  (cell_phi_tag(i,j,k)+cell_phi_tag(i,j,k+1) &
            -cell_phi_tag(i,j-1,k)-cell_phi_tag(i,j-1,k+1)+ &
             cell_phi_tag(i+1,j,k)+cell_phi_tag(i+1,j,k+1) &
            -cell_phi_tag(i+1,j-1,k)-cell_phi_tag(i+1,j-1,k+1))*0.25*dyi


      m3 =  (cell_phi_tag(i,j,k-1)+cell_phi_tag(i,j,k) &
            -cell_phi_tag(i,j-1,k-1)-cell_phi_tag(i,j-1,k)+ &
             cell_phi_tag(i+1,j,k-1)+cell_phi_tag(i+1,j,k) &
            -cell_phi_tag(i+1,j-1,k-1)-cell_phi_tag(i+1,j-1,k))*0.25*dyi

      m4 =  (cell_phi_tag(i,j+1,k-1)+cell_phi_tag(i,j+1,k) &
            -cell_phi_tag(i,j,k-1)-cell_phi_tag(i,j,k)+ &
             cell_phi_tag(i+1,j+1,k-1)+cell_phi_tag(i+1,j+1,k) &
            -cell_phi_tag(i+1,j,k-1)-cell_phi_tag(i+1,j,k))*0.25*dyi



      m5 =  (cell_phi_tag(i,j+1,k)+cell_phi_tag(i,j+1,k+1) &
            -cell_phi_tag(i,j,k)-cell_phi_tag(i,j,k+1)+ &
             cell_phi_tag(i-1,j+1,k)+cell_phi_tag(i-1,j+1,k+1) &
            -cell_phi_tag(i-1,j,k)-cell_phi_tag(i-1,j,k+1))*0.25*dyi
             

      m6 =  (cell_phi_tag(i,j,k)+cell_phi_tag(i,j,k+1) &
            -cell_phi_tag(i,j-1,k)-cell_phi_tag(i,j-1,k+1)+ &
             cell_phi_tag(i-1,j,k)+cell_phi_tag(i-1,j,k+1) &
            -cell_phi_tag(i-1,j-1,k)-cell_phi_tag(i-1,j-1,k+1))*0.25*dyi


      m7 =  (cell_phi_tag(i,j,k-1)+cell_phi_tag(i,j,k) &
            -cell_phi_tag(i,j-1,k-1)-cell_phi_tag(i,j-1,k)+ &
             cell_phi_tag(i-1,j,k-1)+cell_phi_tag(i-1,j,k) &
            -cell_phi_tag(i-1,j-1,k-1)-cell_phi_tag(i-1,j-1,k))*0.25*dyi

      m8 =  (cell_phi_tag(i,j+1,k-1)+cell_phi_tag(i,j+1,k) &
            -cell_phi_tag(i,j,k-1)-cell_phi_tag(i,j,k)+ &
             cell_phi_tag(i-1,j+1,k-1)+cell_phi_tag(i-1,j+1,k) &
            -cell_phi_tag(i-1,j,k-1)-cell_phi_tag(i-1,j,k))*0.25*dyi

      m_av_y= 0.125*(m1+m2+m3+m4+m5+m6+m7+m8)

      m1 = 0
      m2 = 0
      m3 = 0
      m4 = 0
      m5 = 0
      m6 = 0
      m7 = 0
      m8 = 0



   
      m1 =  (cell_phi_tag(i,j,k+1)+cell_phi_tag(i,j+1,k+1) &
            -cell_phi_tag(i,j,k)-cell_phi_tag(i,j+1,k)+&
             cell_phi_tag(i+1,j,k+1)+cell_phi_tag(i+1,j+1,k+1) &
            -cell_phi_tag(i+1,j,k)-cell_phi_tag(i+1,j+1,k))*0.25*dzi


      m2 =  (cell_phi_tag(i,j-1,k+1)+cell_phi_tag(i,j,k+1) &
            -cell_phi_tag(i,j-1,k)-cell_phi_tag(i,j,k)+ &
             cell_phi_tag(i+1,j-1,k+1)+cell_phi_tag(i+1,j,k+1) &
            -cell_phi_tag(i+1,j-1,k)-cell_phi_tag(i+1,j,k))*0.25*dzi
   

      m3 =  (cell_phi_tag(i,j-1,k)+cell_phi_tag(i,j,k) &
            -cell_phi_tag(i,j-1,k-1)-cell_phi_tag(i,j,k-1)+ &
             cell_phi_tag(i+1,j-1,k)+cell_phi_tag(i+1,j,k) &
            -cell_phi_tag(i+1,j-1,k-1)-cell_phi_tag(i+1,j,k-1))*0.25*dzi


      m4 =  (cell_phi_tag(i,j,k)+cell_phi_tag(i,j+1,k) &
            -cell_phi_tag(i,j,k-1)-cell_phi_tag(i,j+1,k-1)+ &
             cell_phi_tag(i+1,j,k)+cell_phi_tag(i+1,j+1,k) &
            -cell_phi_tag(i+1,j,k-1)-cell_phi_tag(i+1,j+1,k-1))*0.25*dzi

      m5 =  (cell_phi_tag(i,j,k+1)+cell_phi_tag(i,j+1,k+1) &
            -cell_phi_tag(i,j,k)-cell_phi_tag(i,j+1,k)+&
             cell_phi_tag(i-1,j,k+1)+cell_phi_tag(i-1,j+1,k+1) &
            -cell_phi_tag(i-1,j,k)-cell_phi_tag(i-1,j+1,k))*0.25*dzi


      m6 =  (cell_phi_tag(i,j-1,k+1)+cell_phi_tag(i,j,k+1) &
            -cell_phi_tag(i,j-1,k)-cell_phi_tag(i,j,k)+ &
             cell_phi_tag(i-1,j-1,k+1)+cell_phi_tag(i-1,j,k+1) &
            -cell_phi_tag(i-1,j-1,k)-cell_phi_tag(i-1,j,k))*0.25*dzi


      m7 =  (cell_phi_tag(i,j-1,k)+cell_phi_tag(i,j,k) &
            -cell_phi_tag(i,j-1,k-1)-cell_phi_tag(i,j,k-1)+ &
             cell_phi_tag(i-1,j-1,k)+cell_phi_tag(i-1,j,k) &
            -cell_phi_tag(i-1,j-1,k-1)-cell_phi_tag(i-1,j,k-1))*0.25*dzi


      m8 =  (cell_phi_tag(i,j,k)+cell_phi_tag(i,j+1,k) &
            -cell_phi_tag(i,j,k-1)-cell_phi_tag(i,j+1,k-1)+ &
             cell_phi_tag(i-1,j,k)+cell_phi_tag(i-1,j+1,k) &
            -cell_phi_tag(i-1,j,k-1)-cell_phi_tag(i-1,j+1,k-1))*0.25*dzi

      m_av_z= 0.125*(m1+m2+m3+m4+m5+m6+m7+m8)



      normal_denum = sqrt(m_av_x*m_av_x+m_av_y*m_av_y+m_av_z*m_av_z+eps)

      nx_surf(i,j,k) = m_av_x/normal_denum
      ny_surf(i,j,k) = m_av_y/normal_denum
      nz_surf(i,j,k) = m_av_z/normal_denum
      nabs_surf(i,j,k) = sqrt(nx_surf(i,j,k)*nx_surf(i,j,k)+&
                                ny_surf(i,j,k)*ny_surf(i,j,k) +&
                               nz_surf(i,j,k)*nz_surf(i,j,k))




     endif

    enddo
  enddo
enddo
return

end subroutine normal_vectors


Subroutine intersect

implicit none
integer::i,j,k,l
real:: nx, ny, nz, n_abs
real:: eps = 1e-12
real:: step =1e-5*dz
real:: xxx,yyy,zzz
logical :: inside
logical :: confirmation
integer:: lmax
real::distance_ghost_intersect
real:: x_ghost,y_ghost,z_ghost
distance_ghost_intersect= 1000.
lmax=100*int( 2*sqrt(3.)*dz/step)
!***************************************************************************************
! Part one: cell-centered
!***************************************************************************************
do k=1,int(solid_height_ratio*k1)
 if (myid.eq.0) print*, 'Calculating Intersect Ponts at k = ', k
 do j= 1,jmax
   do i= 1,imax
    if (abs(nabs_surf(i,j,k)).gt.1e-12) then
       confirmation = .false.
       nx =   nx_surf(i,j,k)
       ny =   ny_surf(i,j,k)
       nz =   nz_surf(i,j,k)
       n_abs = nabs_surf(i,j,k)
       x_ghost =  (i+coords(1)*imax)*dx-0.5*dx
       y_ghost =  (j+coords(2)*jmax)*dy-0.5*dy
       z_ghost =  (k)*dz - 0.5*dz
       do l=0,lmax
         yyy = y_ghost+l*(ny/(n_abs+eps))*step
         zzz = z_ghost+l*(nz/(n_abs+eps))*step
         xxx = x_ghost+l*(nx/(n_abs+eps))*step
         call Solid_Surface(xxx,yyy,zzz,inside)
         if (.not.inside) then
           x_intersect(i,j,k) = xxx
           y_intersect(i,j,k) = yyy
           z_intersect(i,j,k) = zzz
           confirmation = .true.
          exit
         endif
      enddo
      if (.not.confirmation) then
       print*, '--------------------------------------------------------------------------'
       print*,'Error in detecting intersect point at  i , j ,k=' &
              ,i,j,k,'at processor ',myid, 'where the normal vector components are ', &
               nx,ny,nz,n_abs
       print*, '--------------------------------------------------------------------------'
      endif

       distance_ghost_intersect =  sqrt((x_intersect(i,j,k)-x_ghost)*(x_intersect(i,j,k)-x_ghost)+ &
                                        (y_intersect(i,j,k)-y_ghost)*(y_intersect(i,j,k)-y_ghost)+ &
                                        (z_intersect(i,j,k)-z_ghost)*(z_intersect(i,j,k)-z_ghost))
       if (distance_ghost_intersect.gt.sqrt(dx*dx+dy*dy+dz*dz)) then
        print*, '--------------------------------------------------------------------------'
        print*, ' Error in detecting intersect point  processor   ', &
        myid,' : check IBM.f90 - distance_ghost_intesect is ',distance_ghost_intersect, &
       'where the normal vector components are ', nx,ny,nz,n_abs
        print*, '--------------------------------------------------------------------------'
       endif

    endif
   enddo
  enddo
enddo



return 
end subroutine intersect

Subroutine mirrorpoints
implicit none
integer::i,j,k,l,m,n
real:: nx, ny, nz, n_abs
real:: eps = 1e-12
real:: step,step_2
real:: xxx,yyy,zzz
real:: cell_start_y,cell_end_y,cell_start_z,cell_end_z
real::distance_ghost_intersect


step   = march_step
step_2 = march_step*0.05
!***************************************************************************************
! Part one: cell-centered
!***************************************************************************************
do k=1,int(solid_height_ratio*k1)
if (myid.eq.0) print*, 'Calculating Mirror Points at k = ', k
 do j= 1,jmax
   do i= 1,imax
    if (abs(nabs_surf(i,j,k)).gt.1e-12) then

       yyy = (j+coords(2)*jmax)*dy-0.5*dy
       zzz = (k)*dz - 0.5*dz
       xxx = (i+coords(1)*imax)*dx-0.5*dx

       nx =   nx_surf(i,j,k)
       ny =   ny_surf(i,j,k)
       nz =   nz_surf(i,j,k)
       n_abs = nabs_surf(i,j,k)
       l = 0
       distance_ghost_intersect = sqrt( (x_intersect(i,j,k)-xxx)*(x_intersect(i,j,k)-xxx)+ &
                                        (y_intersect(i,j,k)-yyy)*(y_intersect(i,j,k)-yyy)+ &
                                        (z_intersect(i,j,k)-zzz)*(z_intersect(i,j,k)-zzz))





       if  (distance_ghost_intersect.gt.(sqrt(dx*dx+dy*dy+dz*dz))) then
       print*, '--------------------------------------------------------------------------'
           print*, ' Error: in mirro point detection at cell-center at processor ', &
            myid,' : check IBM.f90 - distance_ghost_intesect is ',distance_ghost_intersect, &
            'where the normal vector components are ', nx,ny,nz,n_abs
       print*, '--------------------------------------------------------------------------'
       endif

       x_mirror(i,j,k) = x_intersect(i,j,k)+(nx/(n_abs+eps))*step
       y_mirror(i,j,k) = y_intersect(i,j,k)+(ny/(n_abs+eps))*step !change for 3D
       z_mirror(i,j,k) = z_intersect(i,j,k)+(nz/(n_abs+eps))*step !change for 3D
       deltan(i,j,k)   = distance_ghost_intersect


       x_IP1(i,j,k)    = x_mirror(i,j,k)+(nx/(n_abs+eps))*step_2
       y_IP1(i,j,k)    = y_mirror(i,j,k)+(ny/(n_abs+eps))*step_2
       z_IP1(i,j,k)    = z_mirror(i,j,k)+(nz/(n_abs+eps))*step_2
       
       x_IP2(i,j,k)    = x_IP1(i,j,k)+(nx/(n_abs+eps))*step_2
       y_IP2(i,j,k)    = y_IP1(i,j,k)+(ny/(n_abs+eps))*step_2
       z_IP2(i,j,k)    = z_IP1(i,j,k)+(nz/(n_abs+eps))*step_2

    endif
   enddo
  enddo
enddo



return 
end subroutine mirrorpoints


Subroutine mirrorpoints_ijk
implicit none
integer::i,j,k,l,m,n
real:: nx, ny, nz, n_abs
real:: eps = 1e-12
!real:: step =1e-4*dz
real:: xxx,yyy,zzz
real::cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z
real::distance_ghost_intesect

do k=1,int(solid_height_ratio*k1)
if (myid.eq.0) print*, 'Calculating Mirro Points indices at k = ', k
 do j= 1,jmax
   do i= 1,imax
   !***************************************************************************************
   ! Part one: cell-centered
   !***************************************************************************************
    if (abs(nabs_surf(i,j,k)).gt.1e-12) then
    ! Mirror points
       xxx = x_mirror(i,j,k)
       yyy = y_mirror(i,j,k)
       zzz = z_mirror(i,j,k)
       do l= 1, kmax
         do m= -5,j1+5
           do n= -5,i1+5
             cell_start_x = (n-1+coords(1)*imax)*dx
             cell_end_x   = (n+coords(1)*imax)*dx
             cell_start_y = (m-1+coords(2)*jmax)*dy
             cell_end_y   = (m+coords(2)*jmax)*dy
             cell_start_z = (l-1)*dz
             cell_end_z   = (l)*dz
             if ((yyy.ge.cell_start_y).and.(yyy.lt.cell_end_y).and.&
                 (zzz.ge.cell_start_z).and.(zzz.lt.cell_end_z).and.&
                 (xxx.ge.cell_start_x).and.(xxx.lt.cell_end_x)) then
                  i_mirror(i,j,k) = n 
                  j_mirror(i,j,k) = m
                  k_mirror(i,j,k) = l
              exit


             endif
            enddo
          enddo
       enddo
       if ((i_mirror(i,j,k).eq.i).and.(j_mirror(i,j,k).eq.j).and.(k_mirror(i,j,k).eq.k)) then
             print*, 'Error: Ghost and mirror point are the same(0)'
       endif
       if ((i_mirror(i,j,k).eq.-1000).or.(j_mirror(i,j,k).eq.-1000).or.(k_mirror(i,j,k).eq.-1000)) then
       print*, '--------------------------------------------------------------------------'
         print*,'Error: no grid point detected for mirror point at center i= ', &
         i, ' j= ',j,' k= ',k, 'at processor ',myid
       print*, '--------------------------------------------------------------------------'
      endif

    ! First interpolation points
       xxx = x_IP1(i,j,k)
       yyy = y_IP1(i,j,k)
       zzz = z_IP1(i,j,k)
       do l= 1, kmax
         do m= -5,j1+5
           do n= -5,i1+5
             cell_start_x = (n-1+coords(1)*imax)*dx
             cell_end_x   = (n+coords(1)*imax)*dx
             cell_start_y = (m-1+coords(2)*jmax)*dy
             cell_end_y   = (m+coords(2)*jmax)*dy
             cell_start_z = (l-1)*dz
             cell_end_z   = (l)*dz
             if ((yyy.ge.cell_start_y).and.(yyy.lt.cell_end_y).and.&
                 (zzz.ge.cell_start_z).and.(zzz.lt.cell_end_z).and.&
                 (xxx.ge.cell_start_x).and.(xxx.lt.cell_end_x)) then
                  i_IP1(i,j,k) = n
                  j_IP1(i,j,k) = m
                  k_IP1(i,j,k) = l
              exit
             endif
            enddo
          enddo
       enddo


       if ((i_IP1(i,j,k).eq.-1000).or.(j_IP1(i,j,k).eq.-1000).or.(k_IP1(i,j,k).eq.-1000)) then
       print*, '--------------------------------------------------------------------------'
         print*,'Error: no grid point detected for the first interpolation  point at center i= '&
         ,i, ' j= ',j,' k= ',k, 'at processor ',myid
       print*, '--------------------------------------------------------------------------'
      endif

    ! Second interpolation points
       xxx = x_IP2(i,j,k)
       yyy = y_IP2(i,j,k)
       zzz = z_IP2(i,j,k)
       do l= 1, kmax
         do m= -5,j1+5
           do n= -5,i1+5
             cell_start_x = (n-1+coords(1)*imax)*dx
             cell_end_x   = (n+coords(1)*imax)*dx
             cell_start_y = (m-1+coords(2)*jmax)*dy
             cell_end_y   = (m+coords(2)*jmax)*dy
             cell_start_z = (l-1)*dz
             cell_end_z   = (l)*dz
             if ((yyy.ge.cell_start_y).and.(yyy.lt.cell_end_y).and.&
                 (zzz.ge.cell_start_z).and.(zzz.lt.cell_end_z).and.&
                 (xxx.ge.cell_start_x).and.(xxx.lt.cell_end_x)) then
                  i_IP2(i,j,k) = n
                  j_IP2(i,j,k) = m
                  k_IP2(i,j,k) = l
              exit
             endif

            enddo
          enddo
       enddo


       if ((i_IP2(i,j,k).eq.-1000).or.(j_IP2(i,j,k).eq.-1000).or.(k_IP2(i,j,k).eq.-1000)) then
       print*, '--------------------------------------------------------------------------'
         print*,'Error: no grid point detected for the second interpolation  point at center i= ', &
         i, ' j= ',j,' k= ',k, 'at processor ',myid
       print*, '--------------------------------------------------------------------------'
      endif
    endif
  enddo
 enddo
enddo
return
end Subroutine mirrorpoints_ijk


Subroutine InterpolationWeights
implicit none
real:: x1(7),y1(7),z1(7),x2(7),y2(7),z2(7)
real:: h1(7),h2(7)
integer:: i_p_1(7),j_p_1(7),k_p_1(7),i_p_2(7),j_p_2(7),k_p_2(7)
integer::ii,i,j,k,ii1,ii2,jj1,jj2,kk1,kk2
real:: xx1,xx2,yy1,yy2,zz1,zz2
logical::cond11,cond21,cond12,cond22
real::contribution1(7),contribution2(7)


do k=1,int(solid_height_ratio*k1)
if (myid.eq.0) print*, 'Calculating Interpolation Weights at k = ', k
 do j= 1,jmax
   do i= 1,imax
   !***************************************************************************************
   ! Part one: cell-centered
   !***************************************************************************************
     !Initialization
     xx1               = 0.
     xx2               = 0.
     yy1               = 0.
     yy2               = 0.
     zz1               = 0.
     zz2               = 0.
     cond11            = .false.
     cond21            = .false.
     cond12            = .false.
     cond22            = .false.
     contribution1(:)  =  0.
     contribution2(:)  =  0.
     x1(:)             =  0.
     y1(:)             =  0.
     z1(:)             =  0.
     x2(:)             =  0.
     y2(:)             =  0.
     z2(:)             =  0.
     h1(:)             =  0.
     h2(:)             =  0.
     i_p_1(:)          =  0
     j_p_1(:)          =  0
     k_p_1(:)          =  0
     i_p_2(:)          =  0
     j_p_2(:)          =  0
     k_p_2(:)          =  0
     if (abs(nabs_surf(i,j,k)).gt.1e-12) then
        xx1  =  x_IP1(i,j,k)
        yy1  =  y_IP1(i,j,k)
        zz1  =  z_IP1(i,j,k)     
        ii1  =  i_IP1(i,j,k)
        jj1  =  j_IP1(i,j,k)
        kk1  =  k_IP1(i,j,k)
 
        xx2  =  x_IP2(i,j,k)
        yy2  =  y_IP2(i,j,k)
        zz2  =  z_IP2(i,j,k)
        ii2  =  i_IP2(i,j,k)
        jj2  =  j_IP2(i,j,k)
        kk2  =  k_IP2(i,j,k)
 
 
        i_p_1(1) = ii1 
        j_p_1(1) = jj1+1
        k_p_1(1) = kk1
 
        i_p_2(1) = ii2
        j_p_2(1) = jj2+1
        k_p_2(1) = kk2
 
        
        i_p_1(2) = ii1
        j_p_1(2) = jj1
        k_p_1(2) = kk1+1
 
 
        i_p_2(2) = ii2
        j_p_2(2) = jj2
        k_p_2(2) = kk2+1
 
 
        
        i_p_1(3) = ii1
        j_p_1(3) = jj1-1
        k_p_1(3) = kk1
        
 
        i_p_2(3) = ii2
        j_p_2(3) = jj2-1
        k_p_2(3) = kk2
 
 
 
 
        i_p_1(4) = ii1
        j_p_1(4) = jj1
        k_p_1(4) = kk1-1
        
 
 
        i_p_2(4) = ii2
        j_p_2(4) = jj2
        k_p_2(4) = kk2-1
 
 
 
        i_p_1(5) = ii1
        j_p_1(5) = jj1
        k_p_1(5) = kk1
        
        
        i_p_2(5) = ii2
        j_p_2(5) = jj2
        k_p_2(5) = kk2
 
        
        i_p_1(6) = ii1-1
        j_p_1(6) = jj1
        k_p_1(6) = kk1


        i_p_2(6) = ii2-1
        j_p_2(6) = jj2
        k_p_2(6) = kk2



        i_p_1(7) = ii1+1
        j_p_1(7) = jj1
        k_p_1(7) = kk1


        i_p_2(7) = ii2+1
        j_p_2(7) = jj2
        k_p_2(7) = kk2

        
        do ii = 1,7
           cond11 = (Level_set(i_p_1(ii),j_p_1(ii),k_p_1(ii)).eq.1)
           cond21 = (nabs_surf(i_p_1(ii),j_p_1(ii),k_p_1(ii)).gt.1e-12)
           cond12 = (Level_set(i_p_2(ii),j_p_2(ii),k_p_2(ii)).eq.1)
           cond22 = (nabs_surf(i_p_2(ii),j_p_2(ii),k_p_2(ii)).gt.1e-12)
 
           if (cond11) contribution1(ii) = 1.
           if (cond12) contribution2(ii) = 1.


            x1(ii)=(i_p_1(ii)+coords(1)*imax)*dx-0.5*dx
            y1(ii)=(j_p_1(ii)+coords(2)*jmax)*dy-0.5*dy
            z1(ii)=(k_p_1(ii)               )*dz-0.5*dz
            x2(ii)=(i_p_2(ii)+coords(1)*imax)*dx-0.5*dx
            y2(ii)=(j_p_2(ii)+coords(2)*jmax)*dy-0.5*dy
            z2(ii)=(k_p_2(ii)               )*dz-0.5*dz
        enddo
        
        
        do ii = 1,7
         h1(ii) = sqrt( (x1(ii)-xx1)*(x1(ii)-xx1)+(y1(ii)-yy1)*(y1(ii)-yy1)+(z1(ii)-zz1)*(z1(ii)-zz1)) 
         h2(ii) = sqrt( (x2(ii)-xx2)*(x2(ii)-xx2)+(y2(ii)-yy2)*(y2(ii)-yy2)+(z2(ii)-zz2)*(z2(ii)-zz2)) 
        enddo
        
        
        do ii = 1,7
         WP1(i,j,k,ii)  = (1./(h1(ii)*h1(ii)))*contribution1(ii)
         WP2(i,j,k,ii)  = (1./(h2(ii)*h2(ii)))*contribution2(ii)
        enddo


        !-------- Exceptional cases ---------------------
        do ii = 1,7
         if ((h1(ii).lt.(1e-8)).and.(contribution1(ii).eq. 1)) then
             WP1(i,j,k,1)  = 0.
             WP1(i,j,k,2)  = 0.
             WP1(i,j,k,3)  = 0.
             WP1(i,j,k,4)  = 0.
             WP1(i,j,k,5)  = 0.
             WP1(i,j,k,6)  = 0.
             WP1(i,j,k,7)  = 0.
             WP1(i,j,k,ii)  = 1.
         endif
         if ((h2(ii).lt.(1e-8)).and.(contribution2(ii).eq. 1)) then
             WP2(i,j,k,1)  = 0.
             WP2(i,j,k,2)  = 0.
             WP2(i,j,k,3)  = 0.
             WP2(i,j,k,4)  = 0.
             WP2(i,j,k,5)  = 0.
             WP2(i,j,k,6)  = 0.
             WP2(i,j,k,7)  = 0.
             WP2(i,j,k,ii)  = 1.
         endif
       enddo
       !-------------------------------------------
     endif

   enddo
  enddo
enddo

return 
end subroutine InterpolationWeights


subroutine modify_normal_vectors
implicit none
integer::i,j,k
real::sgn_y,sgn_z
logical::corne_condm
real::sin_cos_max

do k=1,kmax
 do j=1,jmax
  do i=1,imax
sin_cos_max=max(sin(Rotation_angle),cos(Rotation_angle))
corne_condm= (ny_surf(i,j,k).gt.sin_cos_max).or.(nz_surf(i,j,k).gt.sin_cos_max).or. &
             (ny_surf(i,j,k).lt.(-1.*sin_cos_max)).or.(nz_surf(i,j,k).lt.(-1.*sin_cos_max))
  if ((abs(nabs_surf(i,j,k)).gt.1e-12).and.(.not.corne_condm))  then
   sgn_y = ny_surf(i,j,k)/(abs(ny_surf(i,j,k))+1e-12)
   sgn_z = nz_surf(i,j,k)/(abs(nz_surf(i,j,k))+1e-12)
    if ((sgn_y.lt.1e-12).and.(sgn_z.gt.1e-12)) then
      ny_surf(i,j,k) = sgn_y*sin(Rotation_angle) 
      nz_surf(i,j,k) = sgn_z*cos(Rotation_angle)
      nabs_surf(i,j,k) = 1.0
    elseif ((sgn_y.lt.1e-12).and.(sgn_z.lt.1e-12)) then
      ny_surf(i,j,k) = sgn_y*cos(Rotation_angle)
      nz_surf(i,j,k) = sgn_z*sin(Rotation_angle)   
      nabs_surf(i,j,k) = 1.0
    elseif ((sgn_y.gt.1e-12).and.(sgn_z.lt.1e-12)) then
      ny_surf(i,j,k) = sgn_y*sin(Rotation_angle)
      nz_surf(i,j,k) = sgn_z*cos(Rotation_angle)
      nabs_surf(i,j,k) = 1.0
     elseif ((sgn_y.gt.1e-12).and.(sgn_z.gt.1e-12)) then
      ny_surf(i,j,k) = sgn_y*cos(Rotation_angle)
      nz_surf(i,j,k) = sgn_z*sin(Rotation_angle)
      nabs_surf(i,j,k) = 1.02
    endif
  endif
  enddo
 enddo
enddo


end subroutine modify_normal_vectors





subroutine interpolation_mirror(A,iii,jjj,kkk,B)
implicit none
real,dimension(-5:i1+5,-5:j1+5,-5:k1+5),intent(in)::A
integer,intent(in):: iii,jjj,kkk
real,intent(out):: B
real:: eps = 1e-12
real :: q1,q2,WW1(7),WW2(7),B1(7),B2(7),B_1,B_2
integer:: ii1,jj1,kk1,ii2,jj2,kk2,l
!initialization
q1      = 0.
q2      = 0.
WW1(:)  = 0.
WW2(:)  = 0.
B1(:)   = 0.
B2(:)   = 0.
B_1     = 0.
B_2     = 0.
B       = 0.
  ii1     = i_IP1(iii,jjj,kkk)
  jj1     = j_IP1(iii,jjj,kkk)
  kk1     = k_IP1(iii,jjj,kkk)
  ii2     = i_IP2(iii,jjj,kkk)
  jj2     = j_IP2(iii,jjj,kkk)
  kk2     = k_IP2(iii,jjj,kkk)


  WW1(1) = WP1(iii,jjj,kkk,1)
  WW1(2) = WP1(iii,jjj,kkk,2)
  WW1(3) = WP1(iii,jjj,kkk,3)
  WW1(4) = WP1(iii,jjj,kkk,4)
  WW1(5) = WP1(iii,jjj,kkk,5)
  WW1(6) = WP1(iii,jjj,kkk,6)
  WW1(7) = WP1(iii,jjj,kkk,7)
  WW2(1) = WP2(iii,jjj,kkk,1)
  WW2(2) = WP2(iii,jjj,kkk,2)
  WW2(3) = WP2(iii,jjj,kkk,3)
  WW2(4) = WP2(iii,jjj,kkk,4)
  WW2(5) = WP2(iii,jjj,kkk,5)
  WW2(6) = WP2(iii,jjj,kkk,6)
  WW2(7) = WP2(iii,jjj,kkk,7)
  B1(1)    = A(ii1,jj1+1,kk1)
  B1(2)    = A(ii1,jj1,kk1+1)
  B1(3)    = A(ii1,jj1-1,kk1)
  B1(4)    = A(ii1,jj1,kk1-1)
  B1(5)    = A(ii1,jj1,kk1)
  B1(6)    = A(ii1-1,jj1,kk1)
  B1(7)    = A(ii1+1,jj1,kk1)
  B2(1)    = A(ii2,jj2+1,kk2)
  B2(2)    = A(ii2,jj2,kk2+1)
  B2(3)    = A(ii2,jj2-1,kk2)
  B2(4)    = A(ii2,jj2,kk2-1)
  B2(5)    = A(ii2,jj2,kk2)
  B2(6)    = A(ii2-1,jj2,kk2)
  B2(7)    = A(ii2+1,jj2,kk2)

   do l = 1,7
     q1 =  WW1(l) + q1 
     q2 =  WW2(l) + q2   
    enddo


do l = 1,7
 B_1 = (1./q1)*(WW1(l)*B1(l))+B_1
 B_2 = (1./q2)*(WW2(l)*B2(l))+B_2
enddo
B = 2.*B_1-B_2

end subroutine interpolation_mirror

Subroutine interpolation_2D_velocity(UU,VV,WW,iii,jjj,kkk,U_m,V_m,W_m)
implicit none
integer,intent(in) :: iii,jjj,kkk
real,intent(out) ::U_m, V_m, W_m 
real, dimension(-5:i1+5,-5:j1+5,-5:k1+5), intent(in):: UU ,VV ,WW
real:: eps = 1e-12
real :: q1,q2,WW1(7),WW2(7),U1(7),U2(7),V1(7),V2(7),W1(7),W2(7)
real::  U_1,U_2,V_1,V_2,W_1,W_2
integer:: ii1,jj1,kk1,ii2,jj2,kk2,l
!initialization
q1      = 0.
q2      = 0.
WW1(:)  = 0.
WW2(:)  = 0.
U1(:)   = 0.
U2(:)   = 0.
V1(:)   = 0.
V2(:)   = 0.
W1(:)   = 0.
W2(:)   = 0.
U_1     = 0.
U_2     = 0.
V_1     = 0.
V_2     = 0.
W_1     = 0.
W_2     = 0.

  ii1     = i_IP1(iii,jjj,kkk)
  jj1     = j_IP1(iii,jjj,kkk)
  kk1     = k_IP1(iii,jjj,kkk)
  ii2     = i_IP2(iii,jjj,kkk)
  jj2     = j_IP2(iii,jjj,kkk)
  kk2     = k_IP2(iii,jjj,kkk)


  WW1(1) = WP1(iii,jjj,kkk,1)
  WW1(2) = WP1(iii,jjj,kkk,2)
  WW1(3) = WP1(iii,jjj,kkk,3)
  WW1(4) = WP1(iii,jjj,kkk,4)
  WW1(5) = WP1(iii,jjj,kkk,5)
  WW1(6) = WP1(iii,jjj,kkk,6)
  WW1(7) = WP1(iii,jjj,kkk,7)
  WW2(1) = WP2(iii,jjj,kkk,1)
  WW2(2) = WP2(iii,jjj,kkk,2)
  WW2(3) = WP2(iii,jjj,kkk,3)
  WW2(4) = WP2(iii,jjj,kkk,4)
  WW2(5) = WP2(iii,jjj,kkk,5)
  WW2(6) = WP2(iii,jjj,kkk,6)
  WW2(7) = WP2(iii,jjj,kkk,7)

  U1(1)    = UU(ii1,jj1+1,kk1)
  U1(2)    = UU(ii1,jj1,  kk1+1)
  U1(3)    = UU(ii1,jj1-1,kk1)
  U1(4)    = UU(ii1,jj1,  kk1-1)
  U1(5)    = UU(ii1,jj1,  kk1)
  U1(6)    = UU(ii1-1,jj1,  kk1)
  U1(7)    = UU(ii1+1,jj1,  kk1)
  U2(1)    = UU(ii2,jj2+1,kk2)
  U2(2)    = UU(ii2,jj2,  kk2+1)
  U2(3)    = UU(ii2,jj2-1,kk2) 
  U2(4)    = UU(ii2,jj2,  kk2-1)
  U2(5)    = UU(ii2,jj2,  kk2)
  U2(6)    = UU(ii2-1,jj2,  kk2)
  U2(7)    = UU(ii2+1,jj2,  kk2)

  V1(1)    = VV(ii1,jj1+1,kk1)
  V1(2)    = VV(ii1,jj1,  kk1+1)
  V1(3)    = VV(ii1,jj1-1,kk1)
  V1(4)    = VV(ii1,jj1,  kk1-1)
  V1(5)    = VV(ii1,jj1,  kk1)
  V1(6)    = VV(ii1-1,jj1,  kk1)
  V1(7)    = VV(ii1+1,jj1,  kk1)
  V2(1)    = VV(ii2,jj2+1,kk2)
  V2(2)    = VV(ii2,jj2,  kk2+1)
  V2(3)    = VV(ii2,jj2-1,kk2) 
  V2(4)    = VV(ii2,jj2,  kk2-1)
  V2(5)    = VV(ii2,jj2,  kk2)
  V2(6)    = VV(ii2-1,jj2,  kk2)
  V2(7)    = VV(ii2+1,jj2,  kk2)


  W1(1)    = WW(ii1,jj1+1,kk1)
  W1(2)    = WW(ii1,jj1,  kk1+1)
  W1(3)    = WW(ii1,jj1-1,kk1)
  W1(4)    = WW(ii1,jj1,  kk1-1)
  W1(5)    = WW(ii1,jj1,  kk1)
  W1(6)    = WW(ii1-1,jj1,  kk1)
  W1(7)    = WW(ii1+1,jj1,  kk1)
  W2(1)    = WW(ii2,jj2+1,kk2)
  W2(2)    = WW(ii2,jj2,  kk2+1)
  W2(3)    = WW(ii2,jj2-1,kk2)
  W2(4)    = WW(ii2,jj2,  kk2-1)
  W2(5)    = WW(ii2,jj2,  kk2)
  W2(6)    = WW(ii2-1,jj2,  kk2)
  W2(7)    = WW(ii2+1,jj2,  kk2)




   do l = 1,7
     q1 =  WW1(l) + q1 
     q2 =  WW2(l) + q2   
    enddo


do l = 1,7
 U_1 = (1./q1)*(WW1(l)*U1(l))+U_1
 U_2 = (1./q2)*(WW2(l)*U2(l))+U_2 
 V_1 = (1./q1)*(WW1(l)*V1(l))+V_1
 V_2 = (1./q2)*(WW2(l)*V2(l))+V_2
 W_1 = (1./q1)*(WW1(l)*W1(l))+W_1
 W_2 = (1./q2)*(WW2(l)*W2(l))+W_2

enddo
u_m = 2.*U_1-U_2
V_m = 2.*V_1-V_2
W_m = 2.*W_1-W_2

return
end subroutine interpolation_2D_velocity


Subroutine interpolation_2D_dphi(iii,jjj,kkk,dphidx,dphidy,dphidz)
implicit none
integer,intent(in) :: iii,jjj,kkk
real,intent(out):: dphidx,dphidy,dphidz
integer::ii,l
real:: eps = 1e-12
real :: q1,q2,WW1(7),WW2(7)!,B1(5),B2(5),B_1,B_2
integer:: ii1,jj1,kk1,ii2,jj2,kk2
real:: phi_ip1(7),phi_im1(7), phi_jp1(7),phi_jm1(7),phi_kp1(7),phi_km1(7)
real:: phi_ip2(7),phi_im2(7),phi_jp2(7),phi_jm2(7),phi_kp2(7),phi_km2(7)
real:: phi_x1(7),phi_y1(7),phi_z1(7),phi_x2(7),phi_y2(7),phi_z2(7)
real:: dphidx1,dphidy1,dphidz1,dphidx2,dphidy2,dphidz2
!initialization
q1          = 0.
q2          = 0.
WW1(:)      = 0.
WW2(:)      = 0.
phi_ip1(:)  = 0.
phi_im1(:)  = 0.
phi_jp1(:)  = 0.
phi_jm1(:)  = 0.
phi_kp1(:)  = 0.
phi_km1(:)  = 0.
phi_ip2(:)  = 0.
phi_im2(:)  = 0.
phi_jp2(:)  = 0.
phi_jm2(:)  = 0.
phi_kp2(:)  = 0.
phi_km2(:)  = 0.
phi_x1(:)  = 0.
phi_y1(:)  = 0.
phi_z1(:)  = 0.
phi_x2(:)  = 0.
phi_y2(:)  = 0.
phi_z2(:)  = 0.
dphidx1    = 0.
dphidy1    = 0.
dphidz1    = 0.
dphidx2    = 0.
dphidy2    = 0.
dphidz2    = 0.
  ii1     = i_IP1(iii,jjj,kkk)
  jj1     = j_IP1(iii,jjj,kkk)
  kk1     = k_IP1(iii,jjj,kkk)
  ii2     = i_IP2(iii,jjj,kkk)
  jj2     = j_IP2(iii,jjj,kkk)
  kk2     = k_IP2(iii,jjj,kkk)


  WW1(1) = WP1(iii,jjj,kkk,1)
  WW1(2) = WP1(iii,jjj,kkk,2)
  WW1(3) = WP1(iii,jjj,kkk,3)
  WW1(4) = WP1(iii,jjj,kkk,4)
  WW1(5) = WP1(iii,jjj,kkk,5)
  WW1(6) = WP1(iii,jjj,kkk,6)
  WW1(7) = WP1(iii,jjj,kkk,7)
  WW2(1) = WP2(iii,jjj,kkk,1)
  WW2(2) = WP2(iii,jjj,kkk,2)
  WW2(3) = WP2(iii,jjj,kkk,3)
  WW2(4) = WP2(iii,jjj,kkk,4)
  WW2(5) = WP2(iii,jjj,kkk,5)
  WW2(6) = WP2(iii,jjj,kkk,6)
  WW2(7) = WP2(iii,jjj,kkk,7)

  phi_ip1(1) = 0.5*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1+1,jj1+1,kk1))
  phi_im1(1) = 0.5*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1-1,jj1+1,kk1)) 
  phi_jp1(1) = 0.5*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1,jj1+2,kk1))
  phi_jm1(1) = 0.5*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1,jj1,kk1))
  phi_kp1(1) = 0.5*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1,jj1+1,kk1+1))
  phi_km1(1) = 0.5*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1,jj1+1,kk1-1))


  phi_ip1(2) = 0.5*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1+1,jj1,kk1+1))
  phi_im1(2) = 0.5*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1-1,jj1-1,kk1+1))
  phi_jp1(2) = 0.5*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1,jj1+1,kk1+1))
  phi_jm1(2) = 0.5*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1,jj1-1,kk1+1))
  phi_kp1(2) = 0.5*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1,jj1,kk1+2))
  phi_km1(2) = 0.5*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1,jj1,kk1))


  phi_ip1(3) = 0.5*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1+1,jj1-1,kk1))
  phi_im1(3) = 0.5*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1-1,jj1-1,kk1))
  phi_jp1(3) = 0.5*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1,jj1,kk1))
  phi_jm1(3) = 0.5*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1,jj1-2,kk1))
  phi_kp1(3) = 0.5*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1,jj1-1,kk1+1))
  phi_km1(3) = 0.5*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1,jj1-1,kk1-1))

  phi_ip1(4) = 0.5*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1+1,jj1,kk1-1))
  phi_im1(4) = 0.5*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1-1,jj1-1,kk1-1))
  phi_jp1(4) = 0.5*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1,jj1+1,kk1-1))
  phi_jm1(4) = 0.5*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1,jj1-1,kk1-1))
  phi_kp1(4) = 0.5*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1,jj1,kk1))
  phi_km1(4) = 0.5*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1,jj1,kk1-2))


  phi_ip1(5) = 0.5*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1+1,jj1,kk1))
  phi_im1(5) = 0.5*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1-1,jj1,kk1))
  phi_jp1(5) = 0.5*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1,jj1+1,kk1))
  phi_jm1(5) = 0.5*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1,jj1-1,kk1))
  phi_kp1(5) = 0.5*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1,jj1,kk1+1))
  phi_km1(5) = 0.5*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1,jj1,kk1-1))


  phi_ip1(6) = 0.5*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1,jj1,kk1))
  phi_im1(6) = 0.5*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-2,jj1,kk1))
  phi_jp1(6) = 0.5*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-1,jj1+1,kk1))
  phi_jm1(6) = 0.5*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-1,jj1-1,kk1))
  phi_kp1(6) = 0.5*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-1,jj1,kk1+1))
  phi_km1(6) = 0.5*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-1,jj1,kk1-1))


  phi_ip1(7) = 0.5*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+2,jj1,kk1))
  phi_im1(7) = 0.5*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1,jj1,kk1))
  phi_jp1(7) = 0.5*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+1,jj1+1,kk1))
  phi_jm1(7) = 0.5*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+1,jj1-1,kk1))
  phi_kp1(7) = 0.5*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+1,jj1,kk1+1))
  phi_km1(7) = 0.5*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+1,jj1,kk1-1))



  phi_ip2(1) = 0.5*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2+1,jj2+1,kk2))
  phi_im2(1) = 0.5*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2-1,jj2+1,kk2))
  phi_jp2(1) = 0.5*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2,jj2+2,kk2))
  phi_jm2(1) = 0.5*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2,jj2,kk2))
  phi_kp2(1) = 0.5*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2,jj2+1,kk2+1))
  phi_km2(1) = 0.5*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2,jj2+1,kk2-1))


  phi_ip2(2) = 0.5*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2+1,jj2,kk2+1))
  phi_im2(2) = 0.5*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2-1,jj2-1,kk2+1))
  phi_jp2(2) = 0.5*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2,jj2+1,kk2+1))
  phi_jm2(2) = 0.5*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2,jj2-1,kk2+1))
  phi_kp2(2) = 0.5*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2,jj2,kk2+2))
  phi_km2(2) = 0.5*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2,jj2,kk2))


  phi_ip2(3) = 0.5*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2+1,jj2-1,kk2))
  phi_im2(3) = 0.5*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2-1,jj2-1,kk2))
  phi_jp2(3) = 0.5*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2,jj2,kk2))
  phi_jm2(3) = 0.5*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2,jj2-2,kk2))
  phi_kp2(3) = 0.5*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2,jj2-1,kk2+1))
  phi_km2(3) = 0.5*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2,jj2-1,kk2-1))

  phi_ip2(4) = 0.5*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2+1,jj2,kk2-1))
  phi_im2(4) = 0.5*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2-1,jj2-1,kk2-1))
  phi_jp2(4) = 0.5*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2,jj2+1,kk2-1))
  phi_jm2(4) = 0.5*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2,jj2-1,kk2-1))
  phi_kp2(4) = 0.5*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2,jj2,kk2))
  phi_km2(4) = 0.5*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2,jj2,kk2-2))


  phi_ip2(5) = 0.5*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2+1,jj2,kk2))
  phi_im2(5) = 0.5*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2-1,jj2,kk2))
  phi_jp2(5) = 0.5*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2,jj2+1,kk2))
  phi_jm2(5) = 0.5*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2,jj2-1,kk2))
  phi_kp2(5) = 0.5*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2,jj2,kk2+1))
  phi_km2(5) = 0.5*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2,jj2,kk2-1))


  phi_ip2(6) = 0.5*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2,jj2,kk2))
  phi_im2(6) = 0.5*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-2,jj2,kk2))
  phi_jp2(6) = 0.5*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-1,jj2+1,kk2))
  phi_jm2(6) = 0.5*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-1,jj2-1,kk2))
  phi_kp2(6) = 0.5*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-1,jj2,kk2+1))
  phi_km2(6) = 0.5*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-1,jj2,kk2-1))


  phi_ip2(7) = 0.5*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+2,jj2,kk2))
  phi_im2(7) = 0.5*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2,jj2,kk2))
  phi_jp2(7) = 0.5*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+1,jj2+1,kk2))
  phi_jm2(7) = 0.5*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+1,jj2-1,kk2))
  phi_kp2(7) = 0.5*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+1,jj2,kk2+1))
  phi_km2(7) = 0.5*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+1,jj2,kk2-1))



  do ii=1,7
    phi_x1(ii) = (phi_ip1(ii)-phi_im1(ii))*dxi
    phi_y1(ii) = (phi_jp1(ii)-phi_jm1(ii))*dyi
    phi_z1(ii) = (phi_kp1(ii)-phi_km1(ii))*dzi
    phi_x2(ii) = (phi_ip2(ii)-phi_im2(ii))*dxi
    phi_y2(ii) = (phi_jp2(ii)-phi_jm2(ii))*dyi
    phi_z2(ii) = (phi_kp2(ii)-phi_km2(ii))*dzi
    q1 =  WW1(ii) + q1
    q2 =  WW2(ii) + q2
  enddo

do l = 1,7
 dphidx1 = (1./q1)*(WW1(l)*phi_x1(l))+dphidx1
 dphidy1 = (1./q1)*(WW1(l)*phi_y1(l))+dphidy1
 dphidz1 = (1./q1)*(WW1(l)*phi_z1(l))+dphidz1
 dphidx2 = (1./q2)*(WW2(l)*phi_x2(l))+dphidx2
 dphidy2 = (1./q2)*(WW2(l)*phi_y2(l))+dphidy2
 dphidz2 = (1./q2)*(WW2(l)*phi_z2(l))+dphidz2
enddo

dphidx = 2.*dphidx1-dphidx2
dphidy = 2.*dphidy1-dphidy2
dphidz = 2.*dphidz1-dphidz2

return


end subroutine interpolation_2D_dphi


subroutine Penalization_center(uu,vv,ww)
implicit none
real,dimension(-5:i1+5,-5:j1+5,-5:k1+5),intent(inout)::uu,vv,ww
integer::i,j,k
forceytot = 0.
   do k=0,k1
    do j=0,j1
     do i=0,i1
       if (abs(slip_length).gt.1e-12) then
        if (abs(nabs_surf(i,j,k)).lt.1e-12) then
         uu(i,j,k)=cell_phi_tag(i,j,k)*uu(i,j,k)
         vv(i,j,k)=cell_phi_tag(i,j,k)*(vv(i,j,k)-forceytot*dt)
        endif
        if (abs(nabs_surf(i,j,k)).lt.1e-12) then
         ww(i,j,k)=cell_phi_tag(i,j,k)*ww(i,j,k)
        endif
       else
         uu(i,j,k)=cell_phi_tag(i,j,k)*uu(i,j,k)
         vv(i,j,k)=cell_phi_tag(i,j,k)*(vv(i,j,k)-forceytot*dt)
         ww(i,j,k)=cell_phi_tag(i,j,k)*ww(i,j,k)
       endif
      enddo
     enddo
   enddo
end subroutine Penalization_center


subroutine Penalization_face(uu,vv,ww)
implicit none
real,dimension(0:i1,0:j1,0:k1),intent(inout)::uu,vv,ww
integer::i,j,k
forceytot = 0.
   do k=0,k1
    do j=0,j1
     do i=0,i1
       if (slip_length.gt.0) then
          if (.not. ((cell_u_tag(i,j,k).lt.1).and.(cell_u_tag(i,j,k).gt.1e-12))) then
             uu(i,j,k)=cell_u_tag(i,j,k)*uu(i,j,k)
          endif
          if (.not. ((cell_v_tag(i,j,k).lt.1).and.(cell_v_tag(i,j,k).gt.1e-12))) then
             vv(i,j,k)=cell_v_tag(i,j,k)*(vv(i,j,k)-forceytot*dt)
          endif
          if (.not. ((cell_w_tag(i,j,k).lt.1).and.(cell_w_tag(i,j,k).gt.1e-12))) then
              ww(i,j,k)=cell_w_tag(i,j,k)*ww(i,j,k)
          endif
       else
              uu(i,j,k)=cell_u_tag(i,j,k)*uu(i,j,k)
              vv(i,j,k)=cell_v_tag(i,j,k)*(vv(i,j,k)-forceytot*dt)
              ww(i,j,k)=cell_w_tag(i,j,k)*ww(i,j,k)
       
       endif
      enddo
     enddo
   enddo
end subroutine Penalization_face




end module mod_IBM
