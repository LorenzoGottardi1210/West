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
! Yu Jin
!
!-----------------------------------------------------------------------
SUBROUTINE hybrid_kernel_term1234(current_spin, hybrid_kd, sf, iterm)
  !-----------------------------------------------------------------------
  !
  ! iterm == 1: \sum_{v'} (\int v_c \phi_{v'} \phi_{v}) a_{v'}
  ! iterm == 2: \sum_{v'} (\int v_c a_{v'} \phi_{v}) \phi_{v'}
  ! iterm == 3: \sum_{v'} (\int v_c a_{v'} \phi_{v}) a_{v'}
  !             if l_eenac & do_eenac:
  !             \sum_{v'} (\int v_c a_{I,v'} \phi_{v}) a_{J,v'}
  ! iterm == 4: \sum_{v'} (\int v_c a_{v'} a_{v}) \phi_{v'}
  !             if l_eenac & do_eenac:
  !             \sum_{v'} (\int v_c a_{I,v'} a_{J,v}) \phi_{v'}
  !
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : omega
  USE fft_base,              ONLY : dffts
  USE types_coulomb,         ONLY : pot3D
  USE mp,                    ONLY : mp_bcast
  USE fft_at_gamma,          ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,double_invfft_gamma
  USE mp_global,             ONLY : inter_image_comm,my_image_id
  USE pwcom,                 ONLY : npw,npwx,isk,ngk
  USE westcom,               ONLY : nbnd_occ,iuwfc,lrwfc,n_trunc_bands,evc1_all,evc1J_all,l_forces,l_eenac,do_eenac
  USE exx,                   ONLY : exxalfa
  USE buffers,               ONLY : get_buffer
  USE distribution_center,   ONLY : kpt_pool,band_group
  USE wavefunctions,         ONLY : evc,psic
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: current_spin
  LOGICAL, INTENT(IN) :: sf
  INTEGER, INTENT(IN) :: iterm
  COMPLEX(DP), INTENT(INOUT) :: hybrid_kd(npwx,band_group%nlocx)
  !
  ! Workspace
  !
  INTEGER :: lbnd, ibnd, ibndp, jbnd, jbndp, jbnd_end, ir, ig, iks_do
  INTEGER :: current_spin_ikq, ikq, nbndval, flnbndval, nbnd_do
  INTEGER :: dffts_nnr
  COMPLEX(DP), ALLOCATABLE :: aux_hyb(:,:)
  COMPLEX(DP), ALLOCATABLE :: psic2(:), caux(:), gaux(:), gaux2(:), raux(:)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL start_clock('hyb_k1234')
  !
  SELECT CASE(iterm)
  CASE(1,2)
     IF(sf) CALL errore('hybrid_kernel','spin-flip is not supported for term 1 or 2',1)
  CASE(3,4)
  CASE DEFAULT
     CALL errore('hybrid_kernel','invalid term',1)
  END SELECT
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(aux_hyb(npwx,band_group%nloc))
  ALLOCATE(psic2(dffts%nnr))
  ALLOCATE(caux(dffts%nnr))
  ALLOCATE(gaux(npwx))
  ALLOCATE(gaux2(npwx))
  ALLOCATE(raux(dffts%nnr))
  !$acc enter data create(aux_hyb,psic2,caux,gaux,gaux2,raux)
  !
  DO ikq = 1,kpt_pool%nloc
     !
     current_spin_ikq = isk(ikq)
     IF(current_spin_ikq /= current_spin) CYCLE
     !
     IF(sf) THEN
        iks_do = flks(ikq)
     ELSE
        iks_do = ikq
     ENDIF
     !
     nbndval = nbnd_occ(ikq)
     flnbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(ikq)
     !
     ! ... read in GS wavefunctions ikq
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,ikq)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     DO lbnd = 1,nbnd_do,2 ! index to be left
        !
        ibnd = band_group%l2g(lbnd)
        ibndp = ibnd+n_trunc_bands
        !
        !$acc kernels present(raux)
        raux(:) = (0._DP,0._DP)
        !$acc end kernels
        !
        IF(iterm == 3) THEN
           jbnd_end = flnbndval
        ELSE
           jbnd_end = nbndval
        ENDIF
        !
        IF(lbnd < nbnd_do) THEN
           SELECT CASE(iterm)
           CASE(1,2,3)
              CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibndp),evc(:,ibndp+1),psic,'Wave')
           CASE(4)
              IF(l_forces .AND. .NOT. do_eenac) &
              & CALL double_invfft_gamma(dffts,npw,npwx,evc1_all(:,ibnd,iks_do),evc1_all(:,ibnd+1,iks_do),psic,'Wave')
              !
              IF(l_eenac .AND. do_eenac) &
              & CALL double_invfft_gamma(dffts,npw,npwx,evc1J_all(:,ibnd,iks_do),evc1J_all(:,ibnd+1,iks_do),psic,'Wave')
           END SELECT
        ELSE
           SELECT CASE(iterm)
           CASE(1,2,3)
              CALL single_invfft_gamma(dffts,npw,npwx,evc(:,ibndp),psic,'Wave')
           CASE(4)
              IF(l_forces .AND. .NOT. do_eenac) &
              & CALL single_invfft_gamma(dffts,npw,npwx,evc1_all(:,ibnd,iks_do),psic,'Wave')
              !
              IF(l_eenac .AND. do_eenac) &
              & CALL single_invfft_gamma(dffts,npw,npwx,evc1J_all(:,ibnd,iks_do),psic,'Wave')
           END SELECT
        ENDIF
        !
        DO jbnd = 1,jbnd_end-n_trunc_bands ! index to be summed
           !
           jbndp = jbnd+n_trunc_bands
           !
           SELECT CASE(iterm)
           CASE(1)
              CALL double_invfft_gamma(dffts,npw,npwx,evc(:,jbndp),evc1_all(:,jbnd,ikq),psic2,'Wave')
           CASE(2)
              CALL double_invfft_gamma(dffts,npw,npwx,evc1_all(:,jbnd,ikq),evc(:,jbndp),psic2,'Wave')
           CASE(3)
              IF(l_forces .AND. .NOT. do_eenac) &
              & CALL single_invfft_gamma(dffts,npw,npwx,evc1_all(:,jbnd,ikq),psic2,'Wave')
              !
              IF(l_eenac .AND. do_eenac) &
              & CALL double_invfft_gamma(dffts,npw,npwx,evc1_all(:,jbnd,ikq),evc1J_all(:,jbnd,ikq),psic2,'Wave')
           CASE(4)
              CALL double_invfft_gamma(dffts,npw,npwx,evc1_all(:,jbnd,iks_do),evc(:,jbndp),psic2,'Wave')
           END SELECT
           !
           IF(lbnd < nbnd_do) THEN
              !
              !$acc parallel loop present(caux,psic,psic2)
              DO ir = 1,dffts_nnr
                 caux(ir) = psic(ir)*REAL(psic2(ir),KIND=DP)/omega
              ENDDO
              !$acc end parallel
              !
              ! Apply the bare Coulomb potential
              !
              CALL double_fwfft_gamma(dffts,npw,npwx,caux,gaux,gaux2,'Wave')
              !
              !$acc parallel loop present(gaux,pot3D,pot3D%sqvc)
              DO ig = 1,npw
                 gaux(ig) = gaux(ig)*(pot3D%sqvc(ig)**2)
              ENDDO
              !$acc end parallel
              !
              !$acc parallel loop present(gaux2,pot3D,pot3D%sqvc)
              DO ig = 1,npw
                 gaux2(ig) = gaux2(ig)*(pot3D%sqvc(ig)**2)
              ENDDO
              !$acc end parallel
              !
              CALL double_invfft_gamma(dffts,npw,npwx,gaux,gaux2,caux,'Wave')
              !
              IF(iterm == 3) THEN
                 IF(l_forces .AND. .NOT. do_eenac) THEN
                    !$acc parallel loop present(raux,psic2,caux)
                    DO ir = 1,dffts_nnr
                       raux(ir) = raux(ir)+REAL(psic2(ir),KIND=DP)*caux(ir)
                    ENDDO
                    !$acc end parallel
                 ENDIF
                 !
                 IF(l_eenac .AND. do_eenac) THEN
                    !$acc parallel loop present(raux,psic2,caux)
                    DO ir = 1,dffts_nnr
                       raux(ir) = raux(ir)+AIMAG(psic2(ir))*caux(ir)
                    ENDDO
                    !$acc end parallel
                 ENDIF
              ELSE
                 !$acc parallel loop present(raux,psic2,caux)
                 DO ir = 1,dffts_nnr
                    raux(ir) = raux(ir)+AIMAG(psic2(ir))*caux(ir)
                 ENDDO
                 !$acc end parallel
              ENDIF
              !
           ELSE
              !
              !$acc parallel loop present(caux,psic,psic2)
              DO ir = 1,dffts_nnr
                 caux(ir) = REAL(psic(ir),KIND=DP)*REAL(psic2(ir),KIND=DP)/omega
              ENDDO
              !$acc end parallel
              !
              ! Apply the bare Coulomb potential
              !
              CALL single_fwfft_gamma(dffts,npw,npwx,caux,gaux,'Wave')
              !
              !$acc parallel loop present(gaux,pot3D,pot3D%sqvc)
              DO ig = 1,npw
                 gaux(ig) = gaux(ig)*(pot3D%sqvc(ig)**2)
              ENDDO
              !$acc end parallel
              !
              CALL single_invfft_gamma(dffts,npw,npwx,gaux,caux,'Wave')
              !
              IF(iterm == 3) THEN
                 IF(l_forces .AND. .NOT. do_eenac) THEN
                    !$acc parallel loop present(raux,psic2,caux)
                    DO ir = 1,dffts_nnr
                       raux(ir) = raux(ir)+REAL(psic2(ir),KIND=DP)*REAL(caux(ir),KIND=DP)
                    ENDDO
                    !$acc end parallel
                 ENDIF
                 !
                 IF(l_eenac .AND. do_eenac) THEN
                    !$acc parallel loop present(raux,psic2,caux)
                    DO ir = 1,dffts_nnr
                       raux(ir) = raux(ir)+AIMAG(psic2(ir))*REAL(caux(ir),KIND=DP)
                    ENDDO
                    !$acc end parallel
                 ENDIF
              ELSE
                 !$acc parallel loop present(raux,psic2,caux)
                 DO ir = 1,dffts_nnr
                    raux(ir) = raux(ir)+AIMAG(psic2(ir))*REAL(caux(ir),KIND=DP)
                 ENDDO
                 !$acc end parallel
              ENDIF
              !
           ENDIF
           !
        ENDDO
        !
        IF(lbnd < nbnd_do) THEN
           CALL double_fwfft_gamma(dffts,npw,npwx,raux,aux_hyb(:,lbnd),aux_hyb(:,lbnd+1),'Wave')
        ELSE
           CALL single_fwfft_gamma(dffts,npw,npwx,raux,aux_hyb(:,lbnd),'Wave')
        ENDIF
        !
     ENDDO
     !
     IF(iterm == 3) THEN
        !
        ! Recompute nbnd_do for the current spin channel
        !
        nbnd_do = 0
        DO lbnd = 1,band_group%nloc
           ibnd = band_group%l2g(lbnd)+n_trunc_bands
           IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
        ENDDO
        !
     ENDIF
     !
     !$acc parallel loop collapse(2) present(hybrid_kd,aux_hyb)
     DO lbnd = 1,nbnd_do
        DO ig = 1,npw
           hybrid_kd(ig,lbnd) = hybrid_kd(ig,lbnd)-aux_hyb(ig,lbnd)*exxalfa
        ENDDO
     ENDDO
     !$acc end parallel
     !
  ENDDO
  !
  !$acc exit data delete(aux_hyb,psic2,caux,gaux,gaux2,raux)
  DEALLOCATE(aux_hyb)
  DEALLOCATE(psic2)
  DEALLOCATE(caux)
  DEALLOCATE(gaux)
  DEALLOCATE(gaux2)
  DEALLOCATE(raux)
  !
  CALL stop_clock('hyb_k1234')
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE bse_kernel_term4(current_spin, bse_kd4, sf)
  !-----------------------------------------------------------------------
  !
  ! \sum_{v'} (\int W a_{v'} a_{v}) \phi_{v'}
  !
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : omega
  USE fft_base,              ONLY : dffts
  USE noncollin_module,      ONLY : npol
  USE types_coulomb,         ONLY : pot3D_x,pot3D_c
  USE mp,                    ONLY : mp_bcast,mp_sum
  USE fft_at_gamma,          ONLY : single_fwfft_gamma,double_invfft_gamma
  USE mp_global,             ONLY : inter_image_comm,my_image_id,intra_bgrp_comm
  USE pwcom,                 ONLY : npw,npwx,isk,ngk
  USE westcom,               ONLY : ev,dvg,n_pdep_eigen_to_use,nbnd_occ,iuwfc,lrwfc,n_trunc_bands,evc1_all
  USE buffers,               ONLY : get_buffer
  USE distribution_center,   ONLY : kpt_pool,band_group
  USE pdep_db,               ONLY : pdep_db_read
  USE wavefunctions,         ONLY : evc,psic
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: current_spin
  LOGICAL, INTENT(IN) :: sf
  COMPLEX(DP), INTENT(INOUT) :: bse_kd4(npwx,band_group%nlocx)
  !
  ! Workspace
  !
  INTEGER :: lbnd, ibnd, jbnd, jbndp, ir, ig, iks_do, ip
  INTEGER :: current_spin_ikq, ikq, nbndval, nbnd_do
  INTEGER :: dffts_nnr
  REAL(DP) :: factor
  COMPLEX(DP), ALLOCATABLE :: aux_bse4(:,:)
  COMPLEX(DP), ALLOCATABLE :: caux(:), gaux(:), raux(:), tau(:)
  REAL(DP), ALLOCATABLE :: dotp(:)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL start_clock('bse_k4')
  !
  dffts_nnr = dffts%nnr
  !
  CALL pdep_db_read(n_pdep_eigen_to_use,lpara=.FALSE.)
  !$acc enter data copyin(dvg)
  !
  ALLOCATE(aux_bse4(npwx,band_group%nloc))
  ALLOCATE(caux(dffts%nnr))
  ALLOCATE(gaux(npwx))
  ALLOCATE(tau(npwx))
  ALLOCATE(raux(dffts%nnr))
  ALLOCATE(dotp(n_pdep_eigen_to_use))
  !$acc enter data create(aux_bse4,caux,gaux,tau,raux,dotp)
  !
  DO ikq = 1,kpt_pool%nloc
     !
     current_spin_ikq = isk(ikq)
     IF(current_spin_ikq /= current_spin) CYCLE
     !
     IF(sf) THEN
        iks_do = flks(ikq)
     ELSE
        iks_do = ikq
     ENDIF
     !
     nbndval = nbnd_occ(ikq)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(ikq)
     !
     ! ... read in GS wavefunctions ikq
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,ikq)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     DO lbnd = 1,nbnd_do ! index to be left
        !
        ibnd = band_group%l2g(lbnd)
        !
        !$acc kernels present(raux)
        raux(:) = (0._DP,0._DP)
        !$acc end kernels
        !
        DO jbnd = 1,nbndval-n_trunc_bands ! index to be summed
           !
           jbndp = jbnd+n_trunc_bands
           !
           ! product of evc1 and evc1
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc1_all(:,ibnd,iks_do),evc1_all(:,jbnd,iks_do),psic,'Wave')
           !
           !$acc parallel loop present(caux,psic)
           DO ir = 1,dffts_nnr
              caux(ir) = CMPLX(REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))/omega,KIND=DP)
           ENDDO
           !$acc end parallel
           !
           ! Apply the bare Coulomb potential
           !
           CALL single_fwfft_gamma(dffts,npw,npwx,caux,gaux,'Wave')
           !
           !$acc kernels present(tau)
           tau(:) = (0._DP,0._DP)
           !$acc end kernels
           !
           !$acc parallel loop present(tau,gaux,pot3D_x,pot3D_x%sqvc)
           DO ig = 1,npw
              tau(ig) = gaux(ig)*(pot3D_x%sqvc(ig)**2)
           ENDDO
           !$acc end parallel
           !
           !$acc parallel loop present(gaux,pot3D_c,pot3D_c%sqvc)
           DO ig = 1,npw
              gaux(ig) = gaux(ig)*pot3D_c%sqvc(ig)
           ENDDO
           !$acc end parallel
           !
           CALL glbrak_gamma(gaux,dvg,dotp,npw,npwx,1,n_pdep_eigen_to_use,1,npol)
           !
           !$acc update host(dotp)
           !
           CALL mp_sum(dotp,intra_bgrp_comm)
           !
           !$acc kernels present(gaux)
           gaux(:) = (0._DP,0._DP)
           !$acc end kernels
           !
           DO ip = 1,n_pdep_eigen_to_use
              !
              factor = dotp(ip)*ev(ip)/(1._DP-ev(ip))
              !
              !$acc parallel loop present(gaux,dvg)
              DO ig = 1,npw
                 gaux(ig) = gaux(ig)+dvg(ig,ip)*factor
              ENDDO
              !$acc end parallel
              !
           ENDDO
           !
           !$acc parallel loop present(tau,gaux,pot3D_c,pot3D_c%sqvc)
           DO ig = 1,npw
              tau(ig) = tau(ig)+gaux(ig)*pot3D_c%sqvc(ig)
           ENDDO
           !$acc end parallel
           !
           CALL double_invfft_gamma(dffts,npw,npwx,tau,evc(:,jbndp),psic,'Wave')
           !
           !$acc parallel loop present(caux,psic)
           DO ir = 1,dffts_nnr
              caux(ir) = CMPLX(REAL(psic(ir),KIND=DP)*AIMAG(psic(ir)),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           !$acc parallel loop present(raux,caux)
           DO ir = 1,dffts_nnr
              raux(ir) = raux(ir)+caux(ir)
           ENDDO
           !$acc end parallel
           !
        ENDDO
        !
        CALL single_fwfft_gamma(dffts,npw,npwx,raux,aux_bse4(:,lbnd),'Wave')
        !
     ENDDO
     !
     !$acc parallel loop collapse(2) present(bse_kd4,aux_bse4)
     DO lbnd = 1,nbnd_do
        DO ig = 1,npw
           bse_kd4(ig,lbnd) = bse_kd4(ig,lbnd)-aux_bse4(ig,lbnd)
        ENDDO
     ENDDO
     !$acc end parallel
     !
  ENDDO
  !
  !$acc exit data delete(dvg,aux_bse4,caux,gaux,tau,dotp,raux)
  DEALLOCATE(dvg)
  DEALLOCATE(ev)
  DEALLOCATE(aux_bse4)
  DEALLOCATE(caux)
  DEALLOCATE(gaux)
  DEALLOCATE(tau)
  DEALLOCATE(dotp)
  DEALLOCATE(raux)
  !
  CALL stop_clock('bse_k4')
  !
END SUBROUTINE
