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
MODULE plep_io
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    ! ******************************************
    ! WRITE IN G SPACE
    !       wfc is passed distributed in G space
    !       then merged and written in G space
    ! ******************************************
    !
    SUBROUTINE plep_merge_and_write_G(fname,plepg)
      !
      USE kinds,               ONLY : DP,i8b
      USE mp_global,           ONLY : me_bgrp,root_bgrp,nproc_bgrp,intra_bgrp_comm
      USE westcom,             ONLY : nbndval0x,n_trunc_bands
      USE gvect,               ONLY : ig_l2g
      USE pwcom,               ONLY : npwx
      USE noncollin_module,    ONLY : npol
      USE base64_module,       ONLY : islittleendian
      USE west_io,             ONLY : HD_LENGTH,HD_VERSION,HD_ID_VERSION,HD_ID_LITTLE_ENDIAN,HD_ID_DIMENSION
      USE mp_wave,             ONLY : mergewf
      USE mp,                  ONLY : mp_max
      USE distribution_center, ONLY : kpt_pool
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(LEN=*),INTENT(IN) :: fname
      COMPLEX(DP),INTENT(IN) :: plepg(npwx*npol,nbndval0x-n_trunc_bands,kpt_pool%nglob)
      !
      ! Workspae
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: iun,ibnd,ik,npwx_g
      INTEGER :: header(HD_LENGTH)
      INTEGER(i8b) :: offset
      !
      CALL start_clock('plep_write')
      !
      npwx_g = MAXVAL(ig_l2g(1:npwx))
      CALL mp_max(npwx_g,intra_bgrp_comm)
      !
      IF(me_bgrp == root_bgrp) THEN
         !
         header(:) = 0
         header(HD_ID_VERSION) = HD_VERSION
         header(HD_ID_DIMENSION) = npwx_g
         IF(islittleendian()) header(HD_ID_LITTLE_ENDIAN) = 1
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED')
         offset = 1
         WRITE(iun,POS=offset) header
         offset = offset+SIZEOF(header)
         !
      ENDIF
      !
      ! Resume all components
      !
      ALLOCATE(tmp_vec(npwx_g*npol))
      !
      DO ik = 1,kpt_pool%nglob
         DO ibnd = 1,nbndval0x-n_trunc_bands
            !
            tmp_vec(:) = 0._DP
            !
            IF(npol == 2) THEN
               CALL mergewf(plepg(1:npwx,ibnd,ik),tmp_vec(1:npwx_g),npwx,ig_l2g(1:npwx),me_bgrp,&
               & nproc_bgrp,root_bgrp,intra_bgrp_comm)
               CALL mergewf(plepg(npwx+1:npwx*2,ibnd,ik),tmp_vec(npwx_g+1:npwx_g*2),npwx,&
               & ig_l2g(1:npwx),me_bgrp,nproc_bgrp,root_bgrp,intra_bgrp_comm)
            ELSE
               CALL mergewf(plepg(:,ibnd,ik),tmp_vec,npwx,ig_l2g(1:npwx),me_bgrp,nproc_bgrp,&
               & root_bgrp,intra_bgrp_comm)
            ENDIF
            !
            ! ONLY ROOT W/IN BGRP WRITES
            !
            IF(me_bgrp == root_bgrp) THEN
               WRITE(iun,POS=offset) tmp_vec(1:npwx_g*npol)
               offset = offset+SIZEOF(tmp_vec)
            ENDIF
            !
         ENDDO
      ENDDO
      !
      IF(me_bgrp == root_bgrp) CLOSE(iun)
      !
      DEALLOCATE(tmp_vec)
      !
      CALL stop_clock('plep_write')
      !
    END SUBROUTINE
    !
    ! ******************************************
    ! READ IN G SPACE
    !       wfc is read merged in G space
    !       then split in G space
    ! ******************************************
    !
    SUBROUTINE plep_read_G_and_distribute(fname,plepg)
      !
      USE kinds,               ONLY : DP,i8b
      USE mp_global,           ONLY : me_bgrp,root_bgrp,nproc_bgrp,intra_bgrp_comm
      USE westcom,             ONLY : nbndval0x,n_trunc_bands
      USE gvect,               ONLY : ig_l2g
      USE pwcom,               ONLY : npwx
      USE noncollin_module,    ONLY : npol
      USE base64_module,       ONLY : islittleendian
      USE west_io,             ONLY : HD_LENGTH,HD_VERSION,HD_ID_VERSION,HD_ID_LITTLE_ENDIAN,HD_ID_DIMENSION
      USE mp_wave,             ONLY : splitwf
      USE mp,                  ONLY : mp_max
      USE distribution_center, ONLY : kpt_pool
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      CHARACTER(LEN=*),INTENT(IN) :: fname
      COMPLEX(DP),INTENT(OUT) :: plepg(npwx*npol,nbndval0x-n_trunc_bands,kpt_pool%nglob)
      !
      ! Workspace
      !
      COMPLEX(DP),ALLOCATABLE :: tmp_vec(:)
      INTEGER :: iun,ierr,ibnd,ik,npwx_g
      INTEGER :: header(HD_LENGTH)
      INTEGER(i8b) :: offset
      !
      CALL start_clock('plep_read')
      !
      npwx_g = MAXVAL(ig_l2g(1:npwx))
      CALL mp_max(npwx_g,intra_bgrp_comm)
      !
      ! Resume all components
      !
      ALLOCATE(tmp_vec(npwx_g*npol))
      tmp_vec(:) = 0._DP
      plepg(:,:,:) = 0._DP
      !
      IF(me_bgrp == root_bgrp) THEN
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(fname),ACCESS='STREAM',FORM='UNFORMATTED',STATUS='OLD',IOSTAT=ierr)
         IF(ierr /= 0) CALL errore('plep_read','Cannot read file: '//TRIM(fname),1)
         !
         offset = 1
         READ(iun,POS=offset) header
         !
         IF(HD_VERSION /= header(HD_ID_VERSION)) &
         & CALL errore('plep_read','Unknown file format: '//TRIM(fname),1)
         IF(npwx_g /= header(HD_ID_DIMENSION)) &
         & CALL errore('plep_read','Dimension mismatch: '//TRIM(fname),1)
         IF((islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 0)) &
         & .OR. (.NOT. islittleendian() .AND. (header(HD_ID_LITTLE_ENDIAN) == 1))) &
         & CALL errore('plep_read','Endianness mismatch: '//TRIM(fname),1)
         !
         offset = offset+SIZEOF(header)
         !
      ENDIF
      !
      DO ik = 1,kpt_pool%nglob
         DO ibnd = 1,nbndval0x-n_trunc_bands
            !
            ! ONLY ROOT W/IN BGRP READS
            !
            IF(me_bgrp == root_bgrp) THEN
               READ(iun,POS=offset) tmp_vec(1:npwx_g*npol)
               offset = offset+SIZEOF(tmp_vec)
            ENDIF
            !
            IF(npol == 2) THEN
               CALL splitwf(plepg(1:npwx,ibnd,ik),tmp_vec(1:npwx_g),npwx,ig_l2g(1:npwx),me_bgrp,&
               & nproc_bgrp,root_bgrp,intra_bgrp_comm)
               CALL splitwf(plepg(npwx+1:npwx*2,ibnd,ik),tmp_vec(npwx_g+1:npwx_g*2),npwx,&
               & ig_l2g(1:npwx),me_bgrp,nproc_bgrp,root_bgrp,intra_bgrp_comm)
            ELSE
               CALL splitwf(plepg(:,ibnd,ik),tmp_vec,npwx,ig_l2g(1:npwx),me_bgrp,nproc_bgrp,&
               & root_bgrp,intra_bgrp_comm)
            ENDIF
            !
         ENDDO
      ENDDO
      !
      IF(me_bgrp == root_bgrp) CLOSE(iun)
      !
      DEALLOCATE(tmp_vec)
      !
      CALL stop_clock('plep_read')
      !
    END SUBROUTINE
    !
END MODULE
