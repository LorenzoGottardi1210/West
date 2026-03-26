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
SUBROUTINE do_exc_drho()
  !
  USE kinds,                 ONLY : DP
  USE io_push,               ONLY : io_push_title
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE cell_base,             ONLY : omega
  USE fft_base,              ONLY : dffts
  USE pwcom,                 ONLY : npw,npwx,ngk,wg
  USE control_flags,         ONLY : gamma_only
  USE gvect,                 ONLY : gstart
  USE mp,                    ONLY : mp_sum,mp_bcast
  USE mp_global,             ONLY : inter_image_comm,my_image_id,intra_bgrp_comm
  USE buffers,               ONLY : get_buffer
  USE westcom,               ONLY : iuwfc,lrwfc,nbndval0x,nbnd_occ,dvg_exc,evc1_all,westpp_range,&
                                  & westpp_l_spin_flip,westpp_n_liouville_to_use,westpp_save_dir
  USE fft_at_gamma,          ONLY : single_invfft_gamma,double_invfft_gamma
  USE plep_db,               ONLY : plep_db_read
  USE distribution_center,   ONLY : pert,kpt_pool,band_group
  USE class_idistribute,     ONLY : idistribute
  USE types_bz_grid,         ONLY : k_grid
  USE wbse_bgrp,             ONLY : init_gather_bands,gather_bands
  USE west_mp,               ONLY : west_mp_wait
  USE wavefunctions,         ONLY : evc,psic
#if defined(__CUDA)
  USE west_gpu,              ONLY : allocate_gpu,deallocate_gpu
#endif
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: ibnd,jbnd,lbnd,ir,ig,iks,iks_do,iexc,lexc,req,dffts_nnr,nbndval,nbnd_do
  INTEGER :: barra_load
  REAL(DP) :: w1,w2,prod,reduce
  REAL(DP),ALLOCATABLE :: dvgdvg_mat(:,:),aux_r(:),drhox1(:),drhox2(:)
  CHARACTER(LEN=512) :: fname
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  IF(westpp_n_liouville_to_use < 1) CALL errore('do_exc_drho','westpp_n_liouville_to_use < 1',1)
  IF(westpp_range(2) > westpp_n_liouville_to_use) &
  & CALL errore('do_exc_drho','westpp_range(2) > westpp_n_liouville_to_use',1)
  IF(.NOT. gamma_only) &
  & CALL errore('do_exc_drho','unrelaxed differential density requires gamma_only',1)
  !
  ! ... DISTRIBUTE
  !
  pert = idistribute()
  CALL pert%init(westpp_n_liouville_to_use,'i','nvec',.TRUE.)
  kpt_pool = idistribute()
  CALL kpt_pool%init(k_grid%nps,'p','kpt',.FALSE.)
  band_group = idistribute()
  CALL band_group%init(nbndval0x,'b','nbndval',.FALSE.)
  !
  ! READ EIGENVALUES AND VECTORS FROM OUTPUT
  !
  CALL plep_db_read(westpp_n_liouville_to_use)
  !
#if defined(__CUDA)
  CALL allocate_gpu()
#endif
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(dvgdvg_mat(nbndval0x,band_group%nlocx))
  ALLOCATE(drhox1(dffts%nnr))
  ALLOCATE(drhox2(dffts%nnr))
  ALLOCATE(aux_r(dffts%nnr))
  !$acc enter data create(dvgdvg_mat,drhox1,drhox2,aux_r)
  !
  CALL io_push_title('(U)nrelaxed Differential Density')
  !
  CALL init_gather_bands()
  !
  barra_load = 0
  DO lexc = 1,pert%nloc
     iexc = pert%l2g(lexc)
     IF(iexc < westpp_range(1) .OR. iexc > westpp_range(2)) CYCLE
     barra_load = barra_load+1
  ENDDO
  !
  CALL start_bar_type(barra,'westpp',pert%nloc*k_grid%nps)
  !
  DO lexc = 1,pert%nlocx
     !
     ! local -> global
     !
     iexc = pert%l2g(lexc)
     !
     DO iks = 1,kpt_pool%nloc
        !
        !$acc enter data copyin(dvg_exc(:,:,iks,lexc))
        !
        CALL gather_bands(dvg_exc(:,:,iks,lexc),evc1_all(:,:,iks),req)
        !
        IF(westpp_l_spin_flip) THEN
           iks_do = flks(iks)
        ELSE
           iks_do = iks
        ENDIF
        !
        nbndval = nbnd_occ(iks_do)
        !
        nbnd_do = 0
        DO lbnd = 1,band_group%nloc
           ibnd = band_group%l2g(lbnd)
           IF(ibnd > 0 .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
        ENDDO
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
        ! CYCLE here because of mp_bcast above
        !
        IF(iexc < westpp_range(1) .OR. iexc > westpp_range(2)) CYCLE
        !
        ! drhox1
        !
        !$acc kernels present(drhox1)
        drhox1(:) = 0._DP
        !$acc end kernels
        !
        ! double bands @ gamma
        !
        DO lbnd = 1,nbnd_do-MOD(nbnd_do,2),2
           !
           ibnd = band_group%l2g(lbnd)
           jbnd = band_group%l2g(lbnd+1)
           !
           w1 = wg(ibnd,iks_do)/omega
           w2 = wg(jbnd,iks_do)/omega
           !
           CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc(:,lbnd,iks,lexc),&
           & dvg_exc(:,lbnd+1,iks,lexc),psic,'Wave')
           !
           !$acc parallel loop present(drhox1,psic)
           DO ir = 1,dffts_nnr
              drhox1(ir) = drhox1(ir) + w1*REAL(psic(ir),KIND=DP)**2 + w2*AIMAG(psic(ir))**2
           ENDDO
           !$acc end parallel
           !
        ENDDO
        !
        ! single band @ gamma
        !
        IF(MOD(nbnd_do,2) == 1) THEN
           !
           lbnd = nbnd_do
           ibnd = band_group%l2g(lbnd)
           !
           w1 = wg(ibnd,iks_do)/omega
           !
           CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc(:,lbnd,iks,lexc),psic,'Wave')
           !
           !$acc parallel loop present(drhox1,psic)
           DO ir = 1,dffts_nnr
              drhox1(ir) = drhox1(ir) + w1*REAL(psic(ir),KIND=DP)**2
           ENDDO
           !$acc end parallel
           !
        ENDIF
        !
        !$acc update host(drhox1)
        !
        ! < dvg | dvg >
        !
        !$acc kernels present(dvgdvg_mat)
        dvgdvg_mat(:,:) = 0._DP
        !$acc end kernels
        !
        CALL west_mp_wait(req)
#if !defined(__GPU_MPI)
        !$acc update device(evc1_all)
#endif
        !
        !$acc parallel present(dvgdvg_mat,evc1_all,dvg_exc(:,:,iks,lexc))
        !$acc loop collapse(2)
        DO lbnd = 1,nbnd_do
           DO ibnd = 1,nbndval
              !
              reduce = 0._DP
              !$acc loop reduction(+:reduce)
              DO ig = 1,npw
                 reduce = reduce &
                 & + REAL(evc1_all(ig,ibnd,iks),KIND=DP)*REAL(dvg_exc(ig,lbnd,iks,lexc),KIND=DP) &
                 & + AIMAG(evc1_all(ig,ibnd,iks))*AIMAG(dvg_exc(ig,lbnd,iks,lexc))
              ENDDO
              !
              dvgdvg_mat(ibnd,lbnd) = 2._DP*reduce
              !
           ENDDO
        ENDDO
        !$acc end parallel
        !
        IF(gstart == 2) THEN
           !$acc parallel loop collapse(2) present(dvgdvg_mat,evc1_all,dvg_exc(:,:,iks,lexc))
           DO lbnd = 1,nbnd_do
              DO ibnd = 1,nbndval
                 dvgdvg_mat(ibnd,lbnd) = dvgdvg_mat(ibnd,lbnd) &
                 & - REAL(evc1_all(1,ibnd,iks),KIND=DP)*REAL(dvg_exc(1,lbnd,iks,lexc),KIND=DP)
              ENDDO
           ENDDO
           !$acc end parallel
        ENDIF
        !
        !$acc host_data use_device(dvgdvg_mat)
        CALL mp_sum(dvgdvg_mat,intra_bgrp_comm)
        !$acc end host_data
        !
        ! drhox2
        !
        !$acc kernels present(drhox2)
        drhox2(:) = 0._DP
        !$acc end kernels
        !
        DO lbnd = 1,nbnd_do
           !
           ibnd = band_group%l2g(lbnd)
           !
           w1 = wg(ibnd,iks_do)/omega
           !
           CALL single_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),psic,'Wave')
           !
           !$acc parallel loop present(aux_r,psic)
           DO ir = 1,dffts_nnr
              aux_r(ir) = REAL(psic(ir),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           DO jbnd = 1,nbndval,2
              !
              IF(jbnd < nbndval) THEN
                 !
                 CALL double_invfft_gamma(dffts,npw,npwx,evc(:,jbnd),evc(:,jbnd+1),psic,'Wave')
                 !
                 !$acc parallel loop present(aux_r,psic,dvgdvg_mat,drhox2)
                 DO ir = 1,dffts_nnr
                    prod = aux_r(ir) * (REAL(psic(ir),KIND=DP)*dvgdvg_mat(jbnd,lbnd) &
                    &                + AIMAG(psic(ir))*dvgdvg_mat(jbnd+1,lbnd))
                    drhox2(ir) = drhox2(ir) - w1*CMPLX(prod,KIND=DP)
                 ENDDO
                 !$acc end parallel
                 !
              ELSE
                 !
                 CALL single_invfft_gamma(dffts,npw,npwx,evc(:,jbnd),psic,'Wave')
                 !
                 !$acc parallel loop present(aux_r,psic,dvgdvg_mat,drhox2)
                 DO ir = 1,dffts_nnr
                    prod = aux_r(ir) * REAL(psic(ir),KIND=DP)*dvgdvg_mat(jbnd,lbnd)
                    drhox2(ir) = drhox2(ir) - w1*CMPLX(prod,KIND=DP)
                 ENDDO
                 !$acc end parallel
                 !
              ENDIF
              !
           ENDDO
           !
        ENDDO
        !
        !$acc update host(drhox2)
        !
        ! output
        !
        WRITE(fname,'(a,i6.6,a,i6.6)') TRIM(westpp_save_dir)//'/drhoK',iks,'E',iexc
        aux_r(:) = drhox1+drhox2
        CALL dump_r(aux_r,TRIM(fname))
        !
        !$acc exit data delete(dvg_exc(:,:,iks,lexc))
        !
        CALL update_bar_type(barra,'westpp',1)
        !
     ENDDO
     !
  ENDDO
  !
  CALL stop_bar_type(barra,'westpp')
  !
  !$acc exit data delete(dvgdvg_mat,drhox1,drhox2,aux_r)
  DEALLOCATE(dvgdvg_mat)
  DEALLOCATE(drhox1)
  DEALLOCATE(drhox2)
  DEALLOCATE(aux_r)
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
  !
END SUBROUTINE
