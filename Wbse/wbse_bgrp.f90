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
! Victor Yu
!
!-----------------------------------------------------------------------
MODULE wbse_bgrp
  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: init_gather_bands
  PUBLIC :: gather_bands
  !
  INTEGER :: send_count
  INTEGER, ALLOCATABLE :: recv_count(:)
  INTEGER, ALLOCATABLE :: recv_displ(:)
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE init_gather_bands()
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : nbgrp
      USE pwcom,                ONLY : npwx
      USE noncollin_module,     ONLY : npol
      USE westcom,              ONLY : nbndval0x,n_trunc_bands,evc1_all,evc1J_all,l_eenac
      USE distribution_center,  ONLY : kpt_pool,band_group
      !
      IMPLICIT NONE
      !
      ! Workspace
      !
      INTEGER :: nbnd_do,nbnd_loc,ibgrp
      !
      ALLOCATE(recv_count(nbgrp))
      ALLOCATE(recv_displ(nbgrp))
      !
      nbnd_do = nbndval0x-n_trunc_bands
      nbnd_loc = nbnd_do/nbgrp
      send_count = npwx*npol*band_group%nloc
      !
      DO ibgrp = 0,nbgrp-1
         IF(ibgrp < MOD(nbnd_do,nbgrp)) THEN
            recv_count(ibgrp+1) = npwx*npol*(nbnd_loc+1)
            recv_displ(ibgrp+1) = npwx*npol*ibgrp*(nbnd_loc+1)
         ELSE
            recv_count(ibgrp+1) = npwx*npol*nbnd_loc
            recv_displ(ibgrp+1) = npwx*npol*(ibgrp*nbnd_loc+MOD(nbnd_do,nbgrp))
         ENDIF
      ENDDO
      !
      ALLOCATE(evc1_all(npwx*npol,nbnd_do,kpt_pool%nloc))
      !$acc enter data create(evc1_all)
      IF(l_eenac) THEN
         ALLOCATE(evc1J_all(npwx*npol,nbnd_do,kpt_pool%nloc))
         !$acc enter data create(evc1J_all)
      ENDIF
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE gather_bands(distributed,gathered,req)
      !------------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE mp_global,            ONLY : inter_bgrp_comm
      USE pwcom,                ONLY : npwx
      USE noncollin_module,     ONLY : npol
      USE westcom,              ONLY : nbndval0x,n_trunc_bands
      USE distribution_center,  ONLY : band_group
      USE west_mp,              ONLY : west_mp_iallgatherv_start
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(IN) :: distributed(npwx*npol,band_group%nlocx)
      COMPLEX(DP), INTENT(OUT) :: gathered(npwx*npol,nbndval0x-n_trunc_bands)
      INTEGER, INTENT(OUT) :: req
      !
      CALL west_mp_iallgatherv_start(distributed,send_count,gathered,recv_count,recv_displ,&
      & inter_bgrp_comm,req)
      !
    END SUBROUTINE
    !
END MODULE
