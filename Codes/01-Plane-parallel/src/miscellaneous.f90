module miscellaneous
contains


! -----------------------------------------------------------------------------------
subroutine crash(message)
character(LEN=*), intent(in) :: message
write(3,*), message
stop
end subroutine crash
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
subroutine followup(message)
use glbl_prmtrs
character(LEN=*), intent(in) :: message
open(4,file=log_file,position="append")
write(4,*), message
close(4)
end subroutine followup
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
! subroutine eddington(i,roche)
! use glbl_prmtrs
! integer, intent(in) :: i ! if i=1 (resp. 2), we look at the star (resp. BH)
! double precision, intent(out) :: roche
! double precision :: q
!
! if (i==1) q=shpe(q_)!eqpar(GM1_)/eqpar(GM2_)
! if (i==2) q=one/shpe(q_)!eqpar(GM2_)/eqpar(GM1_)
! roche = eqpar(a_) * &
!    ( ( 0.49d0*q**(2.d0/3.d0) ) / (0.6d0*q**(2.d0/3.d0)+dlog(one+q**(1.d0/3.d0))) )
!
! !    ! print*, roche/rsolar, eqpar(Rstr_)/rsolar
! ! if (i==1 .and. eqpar(Rstr_)>roche) then
! !    call crash('Error#04 : the star overflows its Roche lobe')
! ! !!$elseif (i==1 .and. eqpar(Rstr_)<roche) then
! ! !!$   write(2,*), 'Star at ', int(100.d0*(eqpar(Rstr_)/roche)), '% of its Roche lobe'
! ! endif
!
! !!$if (i==2 .and. eqpar(Racc_)>roche) then
! ! Useless now since we consider the minimum w/ 80% of the roche lobe size
! !!$   call crash("Problem : the sphere of influence of radius 10*Racc is larger than its Roche lobe")
! !!$elseif (i==2 .and. eqpar(Racc_)<roche) then
! !!$   write(2,*), 'Sphere of influence at ', int(100.d0*(eqpar(Racc_)/roche)), '% of the BH Roche lobe'
! !!$endif
!
! end subroutine eddington
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
! Computes the time ellapsed from start
subroutine chrono(start,start_mess)
use glbl_prmtrs
double precision, intent(in) :: start
character(len=8), intent(in) :: start_mess ! to have the name of the chrono
double precision :: finish
character(len=200) :: duration

call cpu_time(finish)
write(duration,'(F7.1)') finish-start
duration='It took '//trim(duration)//' sec since '//trim(start_mess)
call followup(duration)

end subroutine chrono
! -----------------------------------------------------------------------------------


end module miscellaneous
