!
! Copyright (C) 2015-2025 M. Govoni
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
SUBROUTINE init_pw_arrays(ncalbec)
  !-----------------------------------------------------------------------
  !
  USE control_flags,          ONLY : gamma_only
  USE scf,                    ONLY : v,vltot,vrs,kedtau
  USE io_files,               ONLY : iunres,tmp_dir,restart_dir
  USE xc_lib,                 ONLY : xclib_dft_is
  USE dfunct,                 ONLY : newd
  USE fft_base,               ONLY : dfftp
  USE mp,                     ONLY : mp_bcast,mp_barrier
  USE mp_global,              ONLY : inter_image_comm,my_image_id
  USE mp_world,               ONLY : world_comm
  USE wavefunctions,          ONLY : evc
  USE klist,                  ONLY : nks
  USE becmod,                 ONLY : becp,allocate_bec_type_acc
  USE uspp,                   ONLY : nkb
  USE noncollin_module,       ONLY : npol
  USE buffers,                ONLY : open_buffer,close_buffer,save_buffer
  USE westcom,                ONLY : n_exx_lowrank,iuwfc,lrwfc
  USE gvecs,                  ONLY : doublegrid
  USE pw_restart_new,         ONLY : read_collected_wfc
  USE lsda_mod,               ONLY : nspin
  USE wvfct,                  ONLY : nbnd,npwx
  USE kinds,                  ONLY : i8b
  USE command_line_options,   ONLY : command_line
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: ncalbec
  !
  ! Workspace
  !
  INTEGER :: ik
  LOGICAL :: exst
  LOGICAL :: exst_mem
  LOGICAL :: l_open_buffer
  LOGICAL :: is_wbse_init
  INTEGER(i8b) :: lrwfc_int8
  INTEGER(i8b) :: recl_int8
  INTEGER(i8b),PARAMETER :: max_int4 = 2147483647
  LOGICAL,EXTERNAL :: matches
  !
  CALL start_clock('init_pw_ar')
  !
  ! Sum self-consistent part (vr) and local part (vltot) of potential
  !
  CALL hinit0()
  !
  CALL set_vrs(vrs,vltot,v%of_r,kedtau,v%kin_r,dfftp%nnr,nspin,doublegrid)
  !
  CALL allocate_bec_type_acc(nkb,ncalbec,becp)
  !
  CALL newd
  !
  iunres = 88
  !
  ! Stop if lrwfc overflows 4-byte integer
  !
  lrwfc_int8 = 1_i8b*nbnd*npwx*npol
  IF(lrwfc_int8 > max_int4) &
  & CALL errore('init_pw_ar','lrwfc too large, rerun with more R&G parallelization',1)
  !
  lrwfc = nbnd*npwx*npol
  !
  ! Decide if buffers are needed
  !
  l_open_buffer = .FALSE.
  IF(.NOT. gamma_only) l_open_buffer = .TRUE.
  IF(nks > 1) l_open_buffer = .TRUE.
  is_wbse_init = matches('wbse_init.x',command_line)
  IF(.NOT. is_wbse_init) THEN
     IF(xclib_dft_is('hybrid') .AND. n_exx_lowrank < 1) l_open_buffer = .TRUE.
  ENDIF
  !
#if defined(__NVCOMPILER) || defined(__PGI)
  IF(l_open_buffer) THEN
     !
     ! Handle nvfortran compiler bug
     ! Fortran direct access I/O fails if record length overflows 4-byte integer
     !
     recl_int8 = 1_i8b*nbnd*npwx*npol*16
     IF(recl_int8 > max_int4) &
     & CALL errore('init_pw_ar','recl too large, rerun with more R&G parallelization',1)
  ENDIF
#endif
  !
  IF(my_image_id == 0) THEN
     IF(l_open_buffer) THEN
        CALL open_buffer(iuwfc,'wfc',lrwfc,1,exst_mem,exst,tmp_dir)
        DO ik = 1,nks
           CALL read_collected_wfc(restart_dir(),ik,evc)
           CALL save_buffer(evc,lrwfc,iuwfc,ik)
        ENDDO
        CALL close_buffer(iuwfc,'KEEP')
     ELSE
        CALL read_collected_wfc(restart_dir(),1,evc)
     ENDIF
  ENDIF
  !
  CALL mp_bcast(evc,0,inter_image_comm)
  !
  ! Ensure wfc files have been written before exx_go
  !
  CALL mp_barrier(world_comm)
  !
  CALL exx_go()
  !
  IF(my_image_id == 0) THEN
     IF(l_open_buffer) CALL open_buffer(iuwfc,'wfc',lrwfc,1,exst_mem,exst,tmp_dir)
  ENDIF
  !
  CALL stop_clock('init_pw_ar')
  !
END SUBROUTINE
