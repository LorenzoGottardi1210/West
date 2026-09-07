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
! Ngoc Linh Nguyen, Victor Yu
!
!---------------------------------------------------------------------
SUBROUTINE wbse_calc_dens(devc, drho, sf)
  !---------------------------------------------------------------------
  !
  ! This subroutine calculates the response charge density
  ! from linear response orbitals and ground state orbitals.
  !
  USE kinds,                  ONLY : DP
  USE cell_base,              ONLY : omega
  USE control_flags,          ONLY : gamma_only
  USE fft_base,               ONLY : dffts
  USE lsda_mod,               ONLY : lsda
  USE noncollin_module,       ONLY : noncolin,npol,nspin_mag,domag
  USE pwcom,                  ONLY : npw,npwx,igk_k,current_k,current_spin,isk,wg,ngk
  USE mp,                     ONLY : mp_sum,mp_bcast
  USE mp_global,              ONLY : my_image_id,inter_image_comm,inter_pool_comm,inter_bgrp_comm
  USE buffers,                ONLY : get_buffer
  USE westcom,                ONLY : iuwfc,lrwfc,nbnd_occ,n_trunc_bands
  USE fft_at_gamma,           ONLY : double_invfft_gamma
  USE fft_at_k,               ONLY : single_fwfft_k,single_invfft_k
  USE distribution_center,    ONLY : kpt_pool,band_group
  USE wavefunctions,          ONLY : evc,psic,psic_nc
#if defined(__CUDA)
  USE west_gpu,               ONLY : tmp_r,psic2,psic_nc2
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: devc(npwx*npol,band_group%nlocx,kpt_pool%nloc)
  ! the perturbed wfcs in reciprocal space
  COMPLEX(DP), INTENT(OUT) :: drho(dffts%nnr,nspin_mag)
  ! the change in charge density matrix (change in (rho_tot, mx, my, mz) if noncolin)
  LOGICAL, INTENT(IN) :: sf
  !
  ! Workspace
  !
  INTEGER :: ir, ibnd, iks, nbndval, lbnd, dffts_nnr, iks_do
  REAL(DP) :: w1
#if !defined(__CUDA)
  REAL(DP), ALLOCATABLE :: tmp_r(:)
  COMPLEX(DP), ALLOCATABLE :: psic2(:)
  COMPLEX(DP), ALLOCATABLE :: psic_nc2(:,:)
#endif
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL start_clock('calc_dens')
  !
  dffts_nnr = dffts%nnr
  !
  !$acc kernels present(drho)
  drho(:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
#if !defined(__CUDA)
  IF(gamma_only) THEN
     ALLOCATE(tmp_r(dffts%nnr))
  ELSE
     IF(noncolin) THEN
        ALLOCATE(psic_nc2(dffts%nnr,npol))
     ELSE
        ALLOCATE(psic2(dffts%nnr))
     ENDIF
  ENDIF
#endif
  !
  DO iks = 1, kpt_pool%nloc  ! KPOINT-SPIN LOOP
     !
     IF(sf) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     ! ... Set k-point and spin
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     !
     ! ... read GS wavefunctions
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     IF(gamma_only) THEN
        !
        !$acc kernels present(tmp_r)
        tmp_r(:) = 0._DP
        !$acc end kernels
        !
        ! double bands @ gamma
        !
        DO lbnd = 1, band_group%nloc
           !
           ibnd = band_group%l2g(lbnd)+n_trunc_bands
           IF(ibnd < 1 .OR. ibnd > nbndval) CYCLE
           !
           w1 = wg(ibnd,iks_do)/omega
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),devc(:,lbnd,iks),psic,'Wave')
           !
           !$acc parallel loop present(tmp_r,psic)
           DO ir = 1, dffts_nnr
              tmp_r(ir) = tmp_r(ir) + w1*REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))
           ENDDO
           !$acc end parallel
           !
        ENDDO
        !
        !$acc parallel loop present(drho,tmp_r)
        DO ir = 1, dffts_nnr
           drho(ir,current_spin) = CMPLX(tmp_r(ir),KIND=DP)
        ENDDO
        !$acc end parallel
        !
     ELSE
        !
        ! only single bands
        !
        DO lbnd = 1, band_group%nloc
           !
           ibnd = band_group%l2g(lbnd)+n_trunc_bands
           IF(ibnd < 1 .OR. ibnd > nbndval) CYCLE
           !
           w1 = wg(ibnd,iks_do)/omega
           !
           IF(noncolin) THEN
              !
              CALL single_invfft_k(dffts,npw,npwx,evc(1:npwx,ibnd),psic_nc(:,1),'Wave',igk_k(:,current_k))
              CALL single_invfft_k(dffts,npw,npwx,evc(npwx+1:npwx*2,ibnd),psic_nc(:,2),'Wave',igk_k(:,current_k))
              CALL single_invfft_k(dffts,npw,npwx,devc(1:npwx,lbnd,iks),psic_nc2(:,1),'Wave',igk_k(:,current_k))
              CALL single_invfft_k(dffts,npw,npwx,devc(npwx+1:npwx*2,lbnd,iks),psic_nc2(:,2),'Wave',igk_k(:,current_k))
              !
              !$acc parallel loop present(drho,psic_nc,psic_nc2)
              DO ir = 1, dffts_nnr
                 drho(ir,1) = drho(ir,1) + w1 * (CONJG(psic_nc(ir,1))*psic_nc2(ir,1) &
                 &                             + CONJG(psic_nc(ir,2))*psic_nc2(ir,2))
              ENDDO
              !$acc end parallel
              !
              IF(domag) THEN
                 !$acc parallel loop present(drho,psic_nc,psic_nc2)
                 DO ir = 1, dffts_nnr
                    drho(ir,2) = drho(ir,2) + w1 * (CONJG(psic_nc(ir,1))*psic_nc2(ir,2) &
                    &                             + CONJG(psic_nc(ir,2))*psic_nc2(ir,1))
                    drho(ir,3) = drho(ir,3) + w1 * (CONJG(psic_nc(ir,1))*psic_nc2(ir,2) &
                    &                             - CONJG(psic_nc(ir,2))*psic_nc2(ir,1)) * (0._DP,-1._DP)
                    drho(ir,4) = drho(ir,4) + w1 * (CONJG(psic_nc(ir,1))*psic_nc2(ir,1) &
                    &                             - CONJG(psic_nc(ir,2))*psic_nc2(ir,2))
                 ENDDO
                 !$acc end parallel
              ENDIF
              !
           ELSE
              !
              CALL single_invfft_k(dffts,npw,npwx,evc(:,ibnd),psic,'Wave',igk_k(:,current_k))
              CALL single_invfft_k(dffts,npw,npwx,devc(:,lbnd,iks),psic2,'Wave',igk_k(:,current_k))
              !
              !$acc parallel loop present(drho,psic,psic2)
              DO ir = 1, dffts_nnr
                 drho(ir,current_spin) = drho(ir,current_spin) + w1*CONJG(psic(ir))*psic2(ir)
              ENDDO
              !$acc end parallel
              !
           ENDIF
           !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  !$acc update host(drho)
  CALL mp_sum(drho,inter_pool_comm)
  CALL mp_sum(drho,inter_bgrp_comm)
  !
#if !defined(__CUDA)
  IF(ALLOCATED(tmp_r)) DEALLOCATE(tmp_r)
  IF(ALLOCATED(psic2)) DEALLOCATE(psic2)
  IF(ALLOCATED(psic_nc2)) DEALLOCATE(psic_nc2)
#endif
  !
  CALL stop_clock('calc_dens')
  !
END SUBROUTINE
