module glbl_prmtrs
  implicit none
  public

double precision, parameter :: smalldble=1.d-10, bigdble=1.d99
integer         , parameter :: intbigdble=100000000
double precision, parameter :: one=1.d0, dpi=dacos(-1.d0), zero=0.d0, half=0.5d0, two=2.d0
character(len=*), parameter :: undefined = 'undefined'

! CGS
double precision, parameter :: kboltz = 1.3806488D-16
double precision, parameter :: stefan = 5.670373d-5
double precision, parameter :: mnucl = 1.6733D-24
double precision, parameter :: msun = 1.989D+33
double precision, parameter :: rsolar = 7.0D+10
double precision, parameter :: secinday= 8.64d4
double precision, parameter :: secinyear= 3.1557600d7
double precision, parameter :: ggrav  = 6.67384D-8
double precision, parameter :: clight = 2.99792458d10
double precision, parameter :: Lsun = 3.846d23
double precision, parameter :: AU = 1.49597871d13
double precision, parameter :: Jsun = 1.1d47 ! ang. mom. for 25d
double precision, parameter :: sigmaE = 6.652458734d-25 ! = sigma_thomson_scattering_on_free_electrons
double precision, parameter :: kpE = 3.98d-1 ! 3.1d-2! = sigma_thomson / mass_proton.
! For kpE, we used the formula (65) of KUDRITZKI+89 w/ I_HE=2 & N_He/N_H = 0.2

INTEGER, PARAMETER :: nshpe=4, nw=4
INTEGER, PARAMETER :: fll_=1, eta_=fll_+1, q_=eta_+1, bet_=q_+1

! integer, parameter :: nshpe=2
! integer, parameter :: q_=1, fll_=q_+1

! INTEGER, PARAMETER :: neqpar=20
! INTEGER, PARAMETER :: aleph_=1, GM1mod_=aleph_+1, GM2_=GM1mod_+1, Om_=GM2_+1, a_=Om_+1, Rstr_=a_+1, Racc_=Rstr_+1, &
!                       CM_=Racc_+1, gm_=CM_+1, &
!                       GM1_=gm_+1, beta_=GM1_+1, vinf_=beta_+1, &
!                       ratioRadiiIni_=vinf_+1, g0lines_=ratioRadiiIni_+1, &
!                       MdotHL_ = g0lines_+1, rout_=MdotHL_+1, egg_=rout_+1, velIni_=egg_+1, roche2_=velIni_+1, NStoL1_=roche2_+1

INTEGER, PARAMETER :: neqpar=5
INTEGER, PARAMETER :: egg_=1, NStoL1_=egg_+1, CM_=NStoL1_+1, potatL1_=CM_+1, Racc_=potatL1_+1

! usr-defined
! Specify the number of lines in the files where the accelerations are stores
double precision, dimension(69) :: rad_andreas, dist_andreas
COMMON rad_andreas, dist_andreas

! GM1_
! R_ --> Rstr_
! L_ <-- Luminosity
! kbT_ <-- Temperature
! Edd_ <-- Gamma_Eddington
! Mdot_ <-- Mdot_0
! rho0_
! cs0_
! al_
! gmRad_

integer, parameter :: NrWind=200000 ! # of radial steps in the subroutine mass_wind
double precision, dimension(2,NrWind) :: glines1D
double precision, dimension(2,NrWind) :: glines
COMMON glines1D, glines
double precision :: th_min, ph_max, ph_min, courantpar
COMMON th_min, ph_max, ph_min, courantpar
double precision :: rmin
COMMON rmin
double precision :: chrono_0, chrono_1, chrono_2, chrono_3, chrono_4, chrono_5
COMMON chrono_0, chrono_1, chrono_2, chrono_3, chrono_4, chrono_5
character(len=8), parameter :: chrono_0_mess='chrono_0', chrono_1_mess='chrono_1', &
chrono_2_mess='chrono_2', chrono_3_mess='chrono_3', chrono_4_mess='chrono_4', chrono_5_mess='chrono_5'

!V!DOUBLE PRECISION :: normt
DOUBLE PRECISION,DIMENSION(nshpe) :: shpe
DOUBLE PRECISION,DIMENSION(neqpar) :: eqpar
!V!DOUBLE PRECISION,DIMENSION(nw+1) :: normvar
COMMON   eqpar,shpe!,normvar,normt,shpe
INTEGER :: i0max, j0max, kmax, Nth_rmp, Nph_rmp
COMMON   i0max, j0max, kmax, Nth_rmp, Nph_rmp
character(len=500) :: folder
COMMON folder
character(len=400) :: mess
COMMON mess
character(len=400) :: log_file, scale_file, par_file, &
big_picture, surf_pts, arrival_map, out_file, err_file, hist_file, &
mass_wind_1D
character(len=400) :: theta_logangle, param_file, info_file, draw_all_in_plane, &
draw_good_in_plane, adim_file, dim_file, parameter_file, file_radiative, file_pressure, Rout, Lorb
character(len=4)   :: law
COMMON law
COMMON log_file, scale_file, par_file, &
big_picture, surf_pts, arrival_map, out_file, err_file, hist_file, &
mass_wind_1D, adim_file, dim_file, file_radiative, file_pressure, Rout, Lorb
COMMON theta_logangle, param_file, info_file, draw_all_in_plane, draw_good_in_plane, parameter_file
logical :: want_it_all, wind_boosted50, fast_wind
COMMON want_it_all, wind_boosted50, fast_wind

end module glbl_prmtrs
