!
! Copyright (C) 2015-2026 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Marco Govoni
!
!-----------------------------------------------------------------------
SUBROUTINE exx_ungo()
  !-----------------------------------------------------------------------
  !
  USE exx,                    ONLY : deallocate_exx
  USE xc_lib,                 ONLY : xclib_dft_is,stop_exx
  USE westcom,                ONLY : westpp_l_compute_tdm,westpp_l_spin_flip,&
                                   & westpp_l_dipole_realspace,code
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  LOGICAL :: do_stopexx
  !
  IF(TRIM(code) == 'WBSE_INIT') THEN
     do_stopexx = .FALSE.
  ELSEIF(TRIM(code) == 'WESTPP') THEN
     do_stopexx = .FALSE.
     IF(westpp_l_compute_tdm .AND. (.NOT. westpp_l_spin_flip) &
     & .AND. (.NOT. westpp_l_dipole_realspace)) do_stopexx = .TRUE.
  ELSE
     do_stopexx = .TRUE.
  ENDIF
  !
  IF(xclib_dft_is('hybrid') .AND. do_stopexx) THEN
     CALL stop_exx()
     CALL deallocate_exx()
  ENDIF
  !
END SUBROUTINE
