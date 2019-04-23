! Compilation - - - - - - -
! gfortran glbl_prmtrs.f90 miscellaneous.f90 mod_kracken.f90 IO.f90 mc.f90 -o mc

program mc
! use init
! use core_module
use IO
use glbl_prmtrs
use miscellaneous

call cpu_time(chrono_0)

!> to index output files (eg log file where the followup subroutine prints)
call give_filenames_v2
!
! call read_arguments()
! call read_par_files()

! print*, 'Initialization starts - - - - - - -'
!
! call initialization
!
! print*, 'Initialization over - - - - - - - -'

print*, 'Core mod starts - - - - - - - - - -'

! call binary_organ_3

print*, 'Terminated - - - - - - - - - - - -'

print*, 'Now, you can run dimensionizing or pipe the adimensioned boundary conditions into VAC'

call chrono(chrono_0,chrono_0_mess)

end program mc
