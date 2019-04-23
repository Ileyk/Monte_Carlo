module IO
  implicit none
  public

contains

! -----------------------------------------------------------------------------------
!> Read the command line arguments passed to amrvac
subroutine read_arguments()
  use mod_kracken
  use glbl_prmtrs

  integer                          :: len, ier, n
  integer                          :: ibegin
  integer                          :: iterm

  ! Specify the options and their default values
  call kracken('cmd','-i amrvac.par -if ' // undefined // &
       ' -slice 0 -collapse 0 --help .false. -convert .false.')

  ! Get the par file(s)
  call retrev('cmd_i', parameter_file, len, ier)

  ! Split the input files, in case multiple were given
  call get_fields_string(parameter_file, " ,'"""//char(9))

end subroutine read_arguments
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
!> Read in the user-supplied parameter-file
subroutine read_par_files()
  use miscellaneous
  ! use draw ! to determine position of L1
  use glbl_prmtrs

  integer :: unitpar=9
  logical            :: file_exists
  ! double precision   :: M1, Rstr, Om, M2
  double precision   :: f, eta, q, beta
  integer            :: Nr

  namelist /andreas/ file_radiative, file_pressure, wind_boosted50

  ! namelist /physical/ M1, Rstr, Om, M2
  namelist /physical/ f, eta, q, law, beta

  namelist /launching/ i0max, th_min, ph_min, ph_max

  namelist /integration/ courantpar

  namelist /mesh/ Rout, Lorb, Nth_rmp, Nph_rmp, Nr

  open(3,file=err_file)

  wind_boosted50 = .false.

  ! Default values
  ! M1   = bigdble
  ! Rstr = bigdble
  ! Om   = bigdble
  ! M2   = bigdble
  f     = bigdble
  eta   = bigdble
  q     = bigdble
  law   = 'noth'
  beta  = bigdble

  i0max  = intbigdble
  th_min = bigdble
  ph_min = bigdble
  ph_max = bigdble

  courantpar = bigdble

  Rout    = ''
  Lorb    = ''
  Nth_rmp = intbigdble
  Nph_rmp = intbigdble
  Nr      = intbigdble

  print *, "Reading " // trim(parameter_file)

  ! Check whether the file exists
  inquire(file=trim(parameter_file), exist=file_exists)

  if (.not. file_exists) &
    call crash("The parameter file " // trim(parameter_file) // " does not exist")

  open(unitpar, file=trim(parameter_file), status='old')
  ! Try to read in the namelists. They can be absent or in a different
  ! order, since we rewind before each read.
  rewind(unitpar)
  read(unitpar, andreas, end=101)

101    rewind(unitpar)
       read(unitpar, physical, end=102)

102    rewind(unitpar)
       read(unitpar, launching, end=103)

103    rewind(unitpar)
       read(unitpar, integration, end=104)

104    rewind(unitpar)
       read(unitpar, mesh, end=105)

105  close(unitpar)

! andreas - - -
! Check whether the files exist
inquire(file=trim(file_radiative), exist=file_exists)
if (.not. file_exists) &
  call crash("The radiative acceleration file " // trim(file_radiative) // " does not exist")
inquire(file=trim(file_pressure), exist=file_exists)
if (.not. file_exists) &
  call crash("The pressure  acceleration file " // trim(file_pressure) // " does not exist")

! physical - - -
if (f>bigdble/two .or. eta>bigdble/two .or. q>bigdble/two .or. law.EQ.'noth') &
  call crash("Missing physical parameters")
if (beta<bigdble/two .and. law .EQ. 'Andr') &
  print*, 'BEWARE, beta exponent provided not used'
shpe(fll_)=f
shpe(eta_)=eta
shpe(q_)=q
shpe(bet_)=beta

! Stellar parameters from Andreas
! eqpar(GM1_)  = M1   ! in solar masses
! eqpar(Rstr_) = Rstr ! in solar radii
! ! Orbital parameters from Andreas
! eqpar(Om_)   = Om   ! period in days
! ! mass of the CO (fiducial)
! ! Because the value of q and f are very sensitive to the value of the orbital separation a
! ! ( = 30 & 87% for 1.75Rstar, = 4.5 & 105% for 1.85Rstar...), we set the mass of the NS
! ! and deduce q, a and f from it. The value we choose is compatible w/ the value of a
! ! within the uncertainties (1.8+or-5%)
! eqpar(GM2_)  = M2   ! in solar masses

! launching - - -
if (i0max>intbigdble/2 .or. th_min>bigdble/two .or. ph_min>bigdble/two .or. ph_max>bigdble/two) &
  call crash("Missing launching parameters")
th_min=dpi*th_min
ph_min=dpi*ph_min
ph_max=dpi*ph_max

if (courantpar>bigdble/two) call crash("Missing integration parameters")

! mesh - - -
! if     (Rout .EQ. "large") then
!   ! Gives L1, the position of the 1st Lagrangian point,
!   ! in units of the norm (ie the orbital sep)
!   ! and in the frame w/ the CM as origin
!   call newton(M1/M2,eqpar(Rout_),1)
!   eqpar(Rout_) = 1.1d0 * dabs(eqpar(Rout_)-(M1/M2)/(one+M1/M2))
! elseif (Rout .EQ. "small") then
!   eqpar(Rout_) = ( ( 0.49d0*(M2/M1)**(2.d0/3.d0) ) / (0.6d0*(M2/M1)**(2.d0/3.d0)+dlog(one+(M2/M1)**(1.d0/3.d0))) )
!   eqpar(Rout_) = 0.9d0 * eqpar(Rout_)
! else
!   call crash("Unknow mode to define the extension of the simulation space")
! endif
! if (Nr>intbigdble/2 .or. Nth_rmp>intbigdble/2 .or. Nph_rmp>intbigdble/2) call crash("Missing mesh parameters")
! if     (Lorb .EQ. 'par' ) then
!   Nth_rmp = Nth_rmp/2 ! since we will work only in the upper part
! elseif (Lorb .EQ. 'prp') then
!   Nph_rmp = Nph_rmp/2 ! since we will work only from pi/2 to 3pi/2
! else
!   call crash("Which orientation for L orbital? par or prp?")
! endif

end subroutine read_par_files
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine give_filenames_v2
use glbl_prmtrs

folder='./output/'

! adim_file=trim(folder)//'parameters'
param_file=trim(folder)//'param' ! where you store the 4 degrees of freedom
info_file=trim(folder)//'info'       ! computed a, q, fll, 1st Andreas radius

! folders
log_file=trim(folder)//'log'
scale_file=trim(folder)//'scale'
! par_file=trim(folder)//'par'
big_picture=trim(folder)//'rk4_traj'
surf_pts=trim(folder)//'pts'
arrival_map=trim(folder)//'map'
! out_file=trim(folder)//'out'
err_file=trim(folder)//'err'
hist_file=trim(folder)//'vel_at_ROCHE_surf_hist'
mass_wind_1D=trim(folder)//'mass_wind_1D'
theta_logangle=trim(folder)//'theta_mesh'

draw_all_in_plane=trim(folder)//'draw_all_in_plane'
draw_good_in_plane=trim(folder)//'draw_good_in_plane'

! resultats=trim(folder)//'results'

! the updated version of out, directuly used in dimensionizing.f90 in vualatv
dim_file =trim(folder)//'dim'  ! where you store the dimensioned computed quantities (but still in (GM2,roche1,L) units)

! Since lines are appended in log_file each time (see in miscellaneous.f90),
! we need to erase the previous one first
call system("rm -f "//log_file)
! et tant qu a faire...
call system("rm -f "//err_file)

end subroutine give_filenames_v2
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
!> Routine to find entries in a string
subroutine get_fields_string(line, delims)
  !> The line from which we want to read
  character(len=*), intent(inout) :: line
  !> A string with delimiters. For example delims = " ,'"""//char(9)
  character(len=*), intent(in)    :: delims

  integer :: ixs_start
  integer :: ixs_end
  integer :: ix, ix_prev

  ix_prev = 0

  ! Find the starting point of the next entry (a non-delimiter value)
  ix = verify(line(ix_prev+1:), delims)

  ixs_start = ix_prev + ix ! This is the absolute position in 'line'

  ! Get the end point of the current entry (next delimiter index minus one)
  ix = scan(line(ixs_start+1:), delims) - 1

  if (ix == -1) then              ! If there is no last delimiter,
    ixs_end = len(line) ! the end of the line is the endpoint
  else
    ixs_end = ixs_start + ix
  end if

  line = line(ixs_start:ixs_end)

end subroutine get_fields_string
! -----------------------------------------------------------------------------------

end module IO
