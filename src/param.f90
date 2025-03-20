module mod_param
use decomp_2d
implicit none
integer, parameter :: ndims = 2
integer, dimension(ndims), parameter :: dims = (/1,32/)
integer,parameter:: number_of_processors=dims(1)*dims(2)
integer, parameter :: itot =4 , jtot = 640 , ktot = 320
integer, parameter :: it1 = itot+1, jt1 = jtot+1, kt1 = ktot+1
integer, parameter :: imax = itot/dims(1), jmax = jtot/dims(2),kmax = ktot
integer, parameter :: i1 = imax+1, j1 = jmax+1, k1 = kmax+1
real, parameter :: lx = itot*(244.0/jtot),ly = 244.0,lz = 122.0
real, parameter :: dxi = itot/lx, dyi = jtot/ly, dzi = ktot/lz
real, parameter :: dx = 1./dxi, dy = 1./dyi, dz = 1./dzi
real, parameter :: pi = acos(-1.)
real, parameter :: picon = acos(-1.)
character(len=3), parameter :: iniu = 'wet' !
logical, parameter :: isfreeslip = .false.
real, parameter :: bulk_v_sup = 1.
character(len=7), parameter :: PhaseField= 'Droplet'! 'Couette'
!***********************************************************************************************
!Non-dimensional parameters
logical, parameter:: Phi01 = .false.
real,parameter:: Reb = 4.0
real,parameter:: Weber = 100.
real,parameter::  Ca= 1.0 !Weber/Reb
real,parameter::  Pec=1.e4
real, parameter:: Cahn_number= 0.022*2.
real,parameter :: Bo= 0.08
!input parameters
real, parameter:: droplet_radius= ly/15.
real,parameter:: PFM_sigma_ref= 0.0107
real,parameter :: u_ref = 1.0
! Calculated parameters
real,parameter::  PFM_l= droplet_radius * Cahn_number
real,parameter:: M_star=0.0005!2*Cahn_number
!real,parameter :: mobility!=2*sqrt(2.)*M_star*u_ref*droplet_radius*droplet_radius/(3*PFM_sigma_ref)
real,parameter :: mobility = 2*sqrt(2.)*u_ref*PFM_l*droplet_radius/(3*Pec*PFM_sigma_ref)
real, parameter:: visc= 3*PFM_sigma_ref*Ca/(2*sqrt(2.)*u_ref)
real,parameter:: PFM_rho= visc*Reb/(u_ref*droplet_radius)
real,parameter:: rho_2= PFM_rho !1.0 !PFM_rho!inside the droplet
real,parameter::vis_2=visc !inside the droplet
real,parameter::  PFM_omega = u_ref/droplet_radius
real,parameter:: relaxtaion = 200. !1./Cahn_number
!**********************
real,parameter:: rho_1=rho_2 *0.01 ! (1.2)/(998)
real,parameter::vis_1=vis_2 *0.01 ! *(1.6e-5)/(1e-3)
real,parameter::phi_bound = -1.
real, parameter:: PFM_phi_range=1.0
!*********************************************************************************************
real,parameter::PFM_thetta= 45.*pi/180
!real,parameter::PFM_mu_f=3*PFM_sigma_ref*droplet_radius/(2*sqrt(2.)*u_ref*relaxtaion*PFM_l)  
real,parameter::PFM_mu_f= 0.!3*PFM_sigma_ref/(2*sqrt(2.)*u_ref*relaxtaion)
real,parameter:: slip_length =  0.25*droplet_radius
real,parameter:: vel_wall_top = 0.!i2*u_ref
real,parameter:: vel_wall_bot = 0.!2*u_ref
logical, parameter :: Wetting = .true.
!**********************************************************************************************
! IBM parameters
character(len=4), parameter :: surface_type = 'sinu' !'Flat' !'RotB'  !'flat' ! 'circ' !'sinu'
!'cosi' 'RotB' 'Crev'
real,parameter:: Rotation_angle =45.0*pi/180.
real,parameter:: Cos_angle = cos(Rotation_angle)
real,parameter:: Sin_angle = sin(Rotation_angle)
real,parameter:: edge_length_of_box =0.5*ly
real,parameter:: h_w = 16. !wall height
real,parameter:: l_w = 50 !wall lentgh
real,parameter:: t_w = 20 !distance between two walls
real,parameter:: h_b = 60
real,parameter:: n_sin = 15-0.5  ! number of sinusodial waves 
real,parameter:: wall_inc = 0   !wall inclination
real,parameter:: march_step = 0.6*sqrt(dy*dy+dz*dz)
real,parameter::solid_height_ratio = 0.75

!***********************************************************************************************
real,parameter:: g_x = 0.0
real,parameter:: g_y = 0.0!(Bo*PFM_sigma_ref/((rho_2-rho_1)*droplet_radius*droplet_radius))*Sin_angle
real,parameter:: g_z = 0.0!(-Bo*PFM_sigma_ref/((rho_2-rho_1)*droplet_radius*droplet_radius))*Cos_angle
!***********************************************************************************************
 character(len=5), parameter :: datadir = 'data/'
 character(len=5), parameter :: partdir = 'part/'
!***********************************************************************************************
end module mod_param
