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
! Stefano Paolo Villani
!
MODULE wbse_nac
!
IMPLICIT NONE
!
CONTAINS
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_nac(dvg_exc_tmp_I, dvg_exc_tmp_J, omega_JI)
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ityp
  USE pwcom,                ONLY : nspin,npwx
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : logfile,nbndval0x,n_trunc_bands,evc1_all,l_genac,l_eenac,do_eenac
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE json_module,          ONLY : json_file
  USE mp_world,             ONLY : mpime,root
  USE io_push,              ONLY : io_push_title
  USE wbse_bgrp,            ONLY : gather_bands
  USE mp,                   ONLY : mp_waitall
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_bse_gpu,deallocate_bse_gpu
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp_I(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(IN), OPTIONAL :: dvg_exc_tmp_J(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  REAL(DP), INTENT(IN), OPTIONAL :: omega_JI
  !
  ! Workspace
  !
  INTEGER :: iks, n, ia, ipol
  INTEGER, ALLOCATABLE :: reqs(:)
  REAL(DP), ALLOCATABLE :: nac_vec(:), dvgdvg_mat(:,:,:), dvgdvg_mat_JI(:,:,:)
  REAL(DP) :: sumnac_vec
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec(:,:,:), zvector(:,:,:), drhox1(:,:), drhox2(:,:)
  TYPE(json_file) :: json
  INTEGER :: iunit
  CHARACTER(LEN=5) :: label
  !
  IF(l_genac .AND. .NOT. do_eenac) THEN
     label = 'genac'
  ELSEIF(l_eenac .AND. do_eenac) THEN
     label = 'eenac'
  ELSE
     CALL errore('wbse_calc_nac','unexpected error',1)
  ENDIF
  !
  CALL start_clock('calc_'//label)
  !
  CALL io_push_title('Compute '//label)
  !
  n = 3 * nat
  !
  ALLOCATE(nac_vec(n))
  nac_vec(:) = 0._DP
  !
  IF(l_eenac .AND. do_eenac) THEN
     !
     ALLOCATE(reqs(kpt_pool%nloc))
     ALLOCATE(dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc))
     ALLOCATE(dvgdvg_mat_JI(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc))
     !$acc enter data create(dvgdvg_mat,dvgdvg_mat_JI)
     ALLOCATE(drhox1(dffts%nnr, nspin))
     !
     DO iks = 1,kpt_pool%nloc
        CALL gather_bands(dvg_exc_tmp_I(:,:,iks), evc1_all(:,:,iks), reqs(iks))
     ENDDO
     !
     ! drhox1
     !
     CALL wbse_calc_drhox1_nac(dvg_exc_tmp_I, dvg_exc_tmp_J, drhox1)
     !
     CALL wbse_nacvec_drhox1(n, dvg_exc_tmp_I, dvg_exc_tmp_J, drhox1, nac_vec, omega_JI)
     !
     ! < dvgI | dvgJ >
     !
     CALL mp_waitall(reqs)
#if !defined(__GPU_MPI)
     !$acc update device(evc1_all)
#endif
     !
     ! For the band parallelization of rhs_zvec_part1:
     ! dvgdvg_mat    (computed between aI and aJ) and
     ! dvgdvg_mat_IJ (computed between aJ and aI)
     !
     CALL wbse_calc_dvgdvg_mat_nac(dvg_exc_tmp_I, dvg_exc_tmp_J, dvgdvg_mat)
     !
     ! To compute dvgdvg_mat_JI:
     ! put the content of dvg_exc_tmp_J into evc1_all
     !
     DO iks = 1,kpt_pool%nloc
        CALL gather_bands(dvg_exc_tmp_J(:,:,iks), evc1_all(:,:,iks), reqs(iks))
     ENDDO
     CALL mp_waitall(reqs)
#if !defined(__GPU_MPI)
     !$acc update device(evc1_all)
#endif
     !
     ! Compute dvgdvg_mat_JI: inputs dvg_exc_tmp_I and dvg_exc_tmp_J switched
     !
     CALL wbse_calc_dvgdvg_mat_nac(dvg_exc_tmp_J, dvg_exc_tmp_I, dvgdvg_mat_JI)
     !
     ! Revert the content of evc1_all back to dvg_exc_tmp_I
     !
     DO iks = 1,kpt_pool%nloc
        CALL gather_bands(dvg_exc_tmp_I(:,:,iks), evc1_all(:,:,iks), reqs(iks))
     ENDDO
     CALL mp_waitall(reqs)
#if !defined(__GPU_MPI)
     !$acc update device(evc1_all)
#endif
     !
     ! drhox2
     !
     ALLOCATE(drhox2(dffts%nnr, nspin))
     !
     CALL wbse_calc_drhox2_nac(dvgdvg_mat, drhox2)
     !
     CALL wbse_nacvec_drhox2(n, dvgdvg_mat, drhox2, nac_vec, omega_JI)
     !
  ENDIF
  !
  ! Z vector
  !
  ALLOCATE(z_rhs_vec(npwx, band_group%nlocx, kpt_pool%nloc))
  ALLOCATE(zvector(npwx, band_group%nlocx, kpt_pool%nloc))
  !$acc enter data create(z_rhs_vec,zvector)
  !
#if defined(__CUDA)
  CALL allocate_bse_gpu(band_group%nlocx)
#endif
  !
  IF(l_genac .AND. .NOT. do_eenac) THEN
     !$acc kernels present(z_rhs_vec,dvg_exc_tmp_I)
     z_rhs_vec(:,:,:) = dvg_exc_tmp_I
     !$acc end kernels
  ENDIF
  !
  IF(l_eenac .AND. do_eenac) THEN
     CALL build_rhs_zvector_eq_eenac(dvg_exc_tmp_I, dvg_exc_tmp_J, dvgdvg_mat, dvgdvg_mat_JI, drhox1, &
     & drhox2, z_rhs_vec, omega_JI)
  ENDIF
  !
  CALL solve_zvector_eq_cg(z_rhs_vec, zvector)
  !
#if defined(__CUDA)
  CALL deallocate_bse_gpu()
#endif
  !
  CALL wbse_nacvec_drhoz_nac(n, zvector, nac_vec)
  !
  !$acc exit data delete(z_rhs_vec,zvector)
  DEALLOCATE(z_rhs_vec)
  DEALLOCATE(zvector)
  !
  CALL io_push_title(label//' total')
  !
  DO ia = 1,nat
     !
     WRITE(stdout, 9035) ia, ityp(ia), label, (nac_vec(3*ia-3+ipol), ipol = 1,3)
     !
  ENDDO
  !
  IF(mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('output.nac_vec.'//label//'_total', nac_vec(1:n))
     !
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     !
     CALL json%destroy()
     !
  ENDIF
  !
  ! enforce total nac_vec to be 0 in each direction
  !
  DO ipol = 1,3
     !
     sumnac_vec = 0._DP
     !
     DO ia = 1,nat
        sumnac_vec = sumnac_vec + nac_vec(3*ia-3+ipol)
     ENDDO
     !
     DO ia = 1,nat
        nac_vec(3*ia-3+ipol) = nac_vec(3*ia-3+ipol) - sumnac_vec/REAL(nat,KIND=DP)
     ENDDO
     !
  ENDDO
  !
  CALL io_push_title(label//' corrected')
  !
  DO ia = 1,nat
     !
     WRITE(stdout, 9035) ia, ityp(ia), label, (nac_vec(3*ia-3+ipol), ipol=1,3)
     !
  ENDDO
  !
  WRITE(stdout,*)
  !
  IF(mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('output.nac_vec.'//label//'_corrected', nac_vec(1:n))
     !
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     !
     CALL json%destroy()
     !
  ENDIF
  !
  DEALLOCATE(nac_vec)
  IF(l_eenac .AND. do_eenac) THEN
     DEALLOCATE(reqs)
     !$acc exit data delete(dvgdvg_mat,dvgdvg_mat_JI)
     DEALLOCATE(dvgdvg_mat)
     DEALLOCATE(dvgdvg_mat_JI)
     DEALLOCATE(drhox1)
     DEALLOCATE(drhox2)
  ENDIF
  !
  CALL stop_clock('calc_'//label)
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   ',A,' = ',3F14.8)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_drhox1_nac(dvg_exc_tmp_I, dvg_exc_tmp_J, drhox1)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE pwcom,                ONLY : isk,lsda,nspin,current_spin,current_k,wg,ngk,npwx,npw
  USE mp,                   ONLY : mp_sum
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : single_invfft_gamma,double_invfft_gamma
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE westcom,              ONLY : nbnd_occ,n_trunc_bands,l_spin_flip
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_pool_comm,inter_bgrp_comm
  USE io_push,              ONLY : io_push_title
  USE wavefunctions,        ONLY : psic
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp_I(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp_J(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(OUT) :: drhox1(dffts%nnr, nspin)
  !
  ! Workspace
  !
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ir, lbnd, ibnd, jbnd, dffts_nnr
  INTEGER :: barra_load
  REAL(DP) :: w1, w2
  REAL(DP), ALLOCATABLE :: tmp_r(:)
  COMPLEX(DP) , ALLOCATABLE :: psic_J(:)
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Compute drhox1')
  !
  dffts_nnr = dffts%nnr
  drhox1(:,:) = (0._DP,0._DP)
  !
  ALLOCATE(tmp_r(dffts%nnr))
  !$acc enter data create(tmp_r)
  !
  ALLOCATE(psic_J(dffts%nnr))
  !$acc enter data create(psic_J)
  !
  barra_load = 0
  !
  DO iks = 1,kpt_pool%nloc
     !
     IF(l_spin_flip) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) barra_load = barra_load+1
     ENDDO
     !
  ENDDO
  !
  CALL start_bar_type(barra,'drhox1',barra_load)
  !
  DO iks = 1,kpt_pool%nloc
     !
     IF(l_spin_flip) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
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
     !$acc kernels present(tmp_r,psic_J)
     tmp_r(:) = 0._DP
     psic_J(:) = 0._DP
     !$acc end kernels
     !
     ! double bands @ gamma
     !
     DO lbnd = 1,nbnd_do-MOD(nbnd_do,2),2
        !
        ibnd = band_group%l2g(lbnd) + n_trunc_bands
        jbnd = band_group%l2g(lbnd+1) + n_trunc_bands
        !
        w1 = wg(ibnd,iks_do)/omega
        w2 = wg(jbnd,iks_do)/omega
        !
        CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp_I(:,lbnd,iks),dvg_exc_tmp_I(:,lbnd+1,iks),psic,'Wave')
        !
        CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp_J(:,lbnd,iks),dvg_exc_tmp_J(:,lbnd+1,iks),psic_J,'Wave')
        !
        !$acc parallel loop present(tmp_r,psic,psic_J)
        DO ir = 1,dffts_nnr
           tmp_r(ir) = tmp_r(ir) + w1*REAL(psic(ir),KIND=DP)*REAL(psic_J(ir),KIND=DP) &
                               & + w2*AIMAG(psic(ir))*AIMAG(psic_J(ir))
        ENDDO
        !$acc end parallel
        !
        CALL update_bar_type(barra,'drhox1',2)
        !
     ENDDO
     !
     ! single band @ gamma
     !
     IF(MOD(nbnd_do,2) == 1) THEN
        !
        lbnd = nbnd_do
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        !
        w1 = wg(ibnd,iks_do)/omega
        !
        CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp_I(:,lbnd,iks),psic,'Wave')
        !
        CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp_J(:,lbnd,iks),psic_J,'Wave')
        !$acc parallel loop present(tmp_r,psic,psic_J)
        DO ir = 1,dffts_nnr
           tmp_r(ir) = tmp_r(ir) + w1*REAL(psic(ir),KIND=DP)*REAL(psic_J(ir),KIND=DP)
        ENDDO
        !$acc end parallel
        !
        CALL update_bar_type(barra,'drhox1',1)
        !
     ENDIF
     !
     !$acc update host(tmp_r)
     !
     drhox1(:,current_spin) = CMPLX(tmp_r,KIND=DP)
     !
  ENDDO
  !
  CALL mp_sum(drhox1,inter_bgrp_comm)
  CALL mp_sum(drhox1,inter_pool_comm)
  !
  CALL stop_bar_type(barra,'drhox1')
  !
  !$acc exit data delete(tmp_r,psic_J)
  DEALLOCATE(tmp_r)
  DEALLOCATE(psic_J)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_nacvec_drhox1(n, dvg_exc_tmp_I, dvg_exc_tmp_J, drhox1, nac_vec, omega_JI)
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ntyp=>nsp,ityp,tau
  USE cell_base,            ONLY : alat,omega
  USE gvect,                ONLY : g,gstart,ngm,ngl,igtongl
  USE uspp,                 ONLY : nkb,vkb
  USE uspp_init,            ONLY : init_us_2
  USE pwcom,                ONLY : isk,igk_k,lsda,nspin,current_spin,current_k,ngk,npwx,npw,xk,wk
  USE mp,                   ONLY : mp_sum
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE westcom,              ONLY : logfile,nbnd_occ,n_trunc_bands,l_spin_flip
  USE vlocal,               ONLY : vloc
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_pool_comm,inter_bgrp_comm,intra_bgrp_comm
  USE json_module,          ONLY : json_file
  USE mp_world,             ONLY : mpime,root
  USE io_push,              ONLY : io_push_title
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_forces_gpu,deallocate_forces_gpu
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: n
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp_I(npwx*npol, band_group%nlocx, kpt_pool%nloc), drhox1(dffts%nnr, nspin)
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp_J(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  REAL(DP), INTENT(IN)  :: omega_JI
  REAL(DP), INTENT(INOUT) :: nac_vec(n)
  !
  ! Workspace
  !
  COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:,:)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, ipol, lbnd, ibnd, ig
  REAL(DP) :: reduce, factor, this_wk
  REAL(DP), ALLOCATABLE :: nacvec_drhox1(:), nacveclc(:,:), rdrhox1(:,:)
  TYPE(json_file) :: json
  INTEGER :: iunit
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Compute nac_vec of drhox1')
  !
  IF(nspin == 2) THEN
     factor = 1._DP
  ELSE
     factor = 0.5_DP
  ENDIF
  !
#if defined(__CUDA)
  CALL allocate_forces_gpu()
#endif
  !
  ALLOCATE(nacvec_drhox1(n))
  ALLOCATE(nacveclc(3, nat))
  ALLOCATE(dvpsi(npwx, band_group%nlocx, 3))
  !$acc enter data create(dvpsi)
  ALLOCATE(rdrhox1(dffts%nnr, nspin))
  !
  nacvec_drhox1(:) = 0._DP
  !
  CALL start_bar_type(barra,'f_drhox1',kpt_pool%nloc*nat)
  !
  ! nonlocal part
  !
  DO iks = 1,kpt_pool%nloc
     !
     IF(l_spin_flip) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
     !
     CALL g2_kin(iks)
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
#if defined(__CUDA)
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.TRUE.)
#else
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.FALSE.)
#endif
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     this_wk = wk(iks)*factor
     !
     DO ia = 1,nat
        !
        ! 1) | dvpsi_i >
        !
        CALL wbse_get_dvpsi_gamma_nonlocal(ia, dvg_exc_tmp_I(:,:,iks), dvpsi)
        !
        ! 2) nacvec_drhox1_i = < dvg | dvpsi_i >
        !
        DO ipol = 1,3
           !
           reduce = 0._DP
           !
           !$acc parallel loop collapse(2) reduction(+:reduce) present(dvg_exc_tmp_J,dvpsi) copy(reduce)
           DO lbnd = 1,nbnd_do
              DO ig = 1,npw
                 reduce = reduce + REAL(dvg_exc_tmp_J(ig,lbnd,iks),KIND=DP)*REAL(dvpsi(ig,lbnd,ipol),KIND=DP) &
                 &               + AIMAG(dvg_exc_tmp_J(ig,lbnd,iks))*AIMAG(dvpsi(ig,lbnd,ipol))
              ENDDO
           ENDDO
           !$acc end parallel
           !
           reduce = 2._DP*reduce
           !
           IF(gstart == 2) THEN
              !$acc parallel loop reduction(+:reduce) present(dvg_exc_tmp_J,dvpsi) copy(reduce)
              DO lbnd = 1,nbnd_do
                 reduce = reduce - REAL(dvg_exc_tmp_J(1,lbnd,iks),KIND=DP)*REAL(dvpsi(1,lbnd,ipol),KIND=DP)
              ENDDO
              !$acc end parallel
           ENDIF
           !
           nacvec_drhox1(3*ia-3+ipol) = nacvec_drhox1(3*ia-3+ipol) + this_wk*reduce
           !
        ENDDO
        !
        CALL update_bar_type(barra,'f_drhox1',1)
        !
     ENDDO
     !
  ENDDO
  !
  CALL mp_sum(nacvec_drhox1,intra_bgrp_comm)
  CALL mp_sum(nacvec_drhox1,inter_bgrp_comm)
  CALL mp_sum(nacvec_drhox1,inter_pool_comm)
  !
  CALL stop_bar_type(barra,'f_drhox1')
  !
  ! local part
  !
  rdrhox1(:,:) = REAL(drhox1,KIND=DP)
  !
  IF(nspin == 2) THEN
     rdrhox1(:,1) = rdrhox1(:,1)+rdrhox1(:,2)
  ENDIF
  !
  CALL force_lc(nat, tau, ityp, ntyp, alat, omega, ngm, ngl, igtongl, g, rdrhox1(:,1), gstart, &
  & gamma_only, vloc, nacveclc)
  !
  nacveclc(:,:) = -factor*nacveclc
  !
  DO ia = 1,nat
     DO ipol = 1,3
        nacvec_drhox1(3*ia-3+ipol) = ( nacvec_drhox1(3*ia-3+ipol) + nacveclc(ipol,ia) ) / (-omega_JI)
     ENDDO
  ENDDO
  !
  nac_vec(:) = nac_vec+nacvec_drhox1
  !
  CALL io_push_title('eenac drhox1')
  !
  DO ia = 1,nat
     !
     WRITE(stdout, 9035) ia, ityp(ia), (nacvec_drhox1(3*ia-3+ipol), ipol = 1,3)
     !
  ENDDO
  !
  WRITE(stdout,*)
  !
  IF(mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('output.nac_vec.eenac_drhox1', nacvec_drhox1(1:n))
     !
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     !
     CALL json%destroy()
     !
  ENDIF
  !
#if defined(__CUDA)
  CALL deallocate_forces_gpu()
#endif
  !
  DEALLOCATE(nacvec_drhox1)
  DEALLOCATE(nacveclc)
  !$acc exit data delete(dvpsi)
  DEALLOCATE(dvpsi)
  DEALLOCATE(rdrhox1)
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   eenac = ',3F14.8)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_dvgdvg_mat_nac(dvg_exc_tmp_I, dvg_exc_tmp_J, dvgdvg_mat)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : gstart
  USE pwcom,                ONLY : ngk,npwx,npw
  USE mp,                   ONLY : mp_sum
  USE noncollin_module,     ONLY : npol
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE westcom,              ONLY : nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip,evc1_all
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : intra_bgrp_comm
  USE io_push,              ONLY : io_push_title
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp_I(npwx*npol, band_group%nlocx, kpt_pool%nloc) !!! It's not really used
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp_J(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  REAL(DP), INTENT(OUT) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc)
  !
  ! Workspace
  !
  INTEGER :: iks, iks_do, nbndval, nbnd_do, lbnd, ibnd, ig
  REAL(DP) :: reduce
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Compute <dvg|dvg>')
  !
  !$acc kernels present(dvgdvg_mat)
  dvgdvg_mat(:,:,:) = 0._DP
  !$acc end kernels
  !
  CALL start_bar_type(barra,'dvgdvg',kpt_pool%nloc)
  !
  DO iks = 1,kpt_pool%nloc
     !
     npw = ngk(iks)
     !
     IF(l_spin_flip) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     !$acc parallel present(dvgdvg_mat,evc1_all,dvg_exc_tmp_J)
     !$acc loop collapse(2)
     DO lbnd = 1,nbnd_do
        DO ibnd = 1,nbndval-n_trunc_bands
           !
           reduce = 0._DP
           !$acc loop reduction(+:reduce)
           DO ig = 1,npw
              ! evc1_all contains dvg_exc_tmp_I
              reduce = reduce &
              & + REAL(evc1_all(ig,ibnd,iks),KIND=DP)*REAL(dvg_exc_tmp_J(ig,lbnd,iks),KIND=DP) &
              & + AIMAG(evc1_all(ig,ibnd,iks))*AIMAG(dvg_exc_tmp_J(ig,lbnd,iks))
           ENDDO
           !
           dvgdvg_mat(ibnd,lbnd,iks) = 2._DP*reduce
           !
        ENDDO
     ENDDO
     !$acc end parallel
     !
     IF(gstart == 2) THEN
        !$acc parallel loop collapse(2) present(dvgdvg_mat,evc1_all,dvg_exc_tmp_J)
        DO lbnd = 1,nbnd_do
           DO ibnd = 1,nbndval - n_trunc_bands
              dvgdvg_mat(ibnd,lbnd,iks) = dvgdvg_mat(ibnd,lbnd,iks) &
              & - REAL(evc1_all(1,ibnd,iks),KIND=DP)*REAL(dvg_exc_tmp_J(1,lbnd,iks),KIND=DP)
           ENDDO
        ENDDO
        !$acc end parallel
     ENDIF
     !
     CALL update_bar_type(barra,'dvgdvg',1)
     !
  ENDDO
  !
  !$acc host_data use_device(dvgdvg_mat)
  CALL mp_sum(dvgdvg_mat,intra_bgrp_comm)
  !$acc end host_data
  !
  CALL stop_bar_type(barra,'dvgdvg')
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_drhox2_nac(dvgdvg_mat, drhox2)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE pwcom,                ONLY : isk,lsda,wg,ngk,current_spin,nspin,npwx,npw
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE fft_base,             ONLY : dffts
  USE fft_at_gamma,         ONLY : double_invfft_gamma,single_invfft_gamma
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,inter_bgrp_comm
  USE io_push,              ONLY : io_push_title
  USE wavefunctions,        ONLY : evc,psic
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(OUT) :: drhox2(dffts%nnr, nspin)
  !
  ! Workspace
  !
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ir, lbnd, ibnd, jbnd, jbndp, dffts_nnr
  INTEGER :: barra_load
  REAL(DP) :: prod, w1
  REAL(DP), ALLOCATABLE :: aux_r(:)
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Compute drhox2')
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(aux_r(dffts%nnr))
  !$acc enter data create(aux_r,drhox2)
  !
  !$acc kernels present(drhox2)
  drhox2(:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
  barra_load = 0
  !
  DO iks = 1,kpt_pool%nloc
     !
     IF(l_spin_flip) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) barra_load = barra_load+1
     ENDDO
     !
  ENDDO
  !
  CALL start_bar_type(barra,'drhox2',barra_load)
  !
  DO iks = 1,kpt_pool%nloc
     !
     IF(l_spin_flip) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Set k-point and spin
     !
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
     DO lbnd = 1,nbnd_do
        !
        ibnd = band_group%l2g(lbnd) + n_trunc_bands
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
        DO jbnd = 1,nbndval-n_trunc_bands,2
           !
           jbndp = jbnd + n_trunc_bands
           !
           IF(jbnd < nbndval-n_trunc_bands) THEN
              !
              CALL double_invfft_gamma(dffts,npw,npwx,evc(:,jbndp),evc(:,jbndp+1),psic,'Wave')
              !
              !$acc parallel loop present(aux_r,psic,dvgdvg_mat,drhox2)
              DO ir = 1,dffts_nnr
                 prod = aux_r(ir) * (REAL(psic(ir),KIND=DP)*dvgdvg_mat(jbnd,lbnd,iks) &
                 &                + AIMAG(psic(ir))*dvgdvg_mat(jbnd+1,lbnd,iks))
                 drhox2(ir,current_spin) = drhox2(ir,current_spin) - w1*CMPLX(prod,KIND=DP)
              ENDDO
              !$acc end parallel
              !
           ELSE
              !
              CALL single_invfft_gamma(dffts,npw,npwx,evc(:,jbndp),psic,'Wave')
              !
              !$acc parallel loop present(aux_r,psic,dvgdvg_mat,drhox2)
              DO ir = 1,dffts_nnr
                 prod = aux_r(ir) * REAL(psic(ir),KIND=DP) * dvgdvg_mat(jbnd,lbnd,iks)
                 drhox2(ir,current_spin) = drhox2(ir,current_spin) - w1*CMPLX(prod,KIND=DP)
              ENDDO
              !$acc end parallel
              !
           ENDIF
           !
        ENDDO
        !
        CALL update_bar_type(barra,'drhox2',1)
        !
     ENDDO
     !
  ENDDO
  !
  !$acc exit data copyout(drhox2)
  !
  CALL mp_sum(drhox2,inter_bgrp_comm)
  CALL mp_sum(drhox2,inter_pool_comm)
  !
  CALL stop_bar_type(barra,'drhox2')
  !
  !$acc exit data delete(aux_r)
  DEALLOCATE(aux_r)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_nacvec_drhox2(n, dvgdvg_mat, drhox2, nac_vec, omega_JI)
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ntyp=>nsp,ityp,tau
  USE cell_base,            ONLY : alat,omega
  USE gvect,                ONLY : g,gstart,ngm,ngl,igtongl
  USE uspp,                 ONLY : nkb,vkb
  USE uspp_init,            ONLY : init_us_2
  USE pwcom,                ONLY : isk,igk_k,lsda,current_spin,nspin,current_k,ngk,npwx,npw,xk,wk
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE fft_base,             ONLY : dffts
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE westcom,              ONLY : iuwfc,lrwfc,logfile,nbnd_occ,nbndval0x,n_trunc_bands,l_spin_flip
  USE vlocal,               ONLY : vloc
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,inter_bgrp_comm,&
                                 & intra_bgrp_comm
  USE json_module,          ONLY : json_file
  USE mp_world,             ONLY : mpime,root
  USE io_push,              ONLY : io_push_title
  USE wavefunctions,        ONLY : evc
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_forces_gpu,deallocate_forces_gpu
  USE cublas
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands, band_group%nlocx, kpt_pool%nloc)
  COMPLEX(DP), INTENT(IN) :: drhox2(dffts%nnr, nspin)
  REAL(DP), INTENT(INOUT) :: nac_vec(n)
  REAL(DP), INTENT(IN)  :: omega_JI
  !
  ! Workspace
  !
  COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:,:), aux1(:,:), aux2(:,:)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, ipol, lbnd, ibnd, ig
  INTEGER :: band_group_myoffset
  REAL(DP) :: reduce, factor, this_wk
  REAL(DP), ALLOCATABLE :: nacvec_drhox2(:), rdrhox2(:,:), nacveclc(:,:)
  TYPE(json_file) :: json
  INTEGER :: iunit
  TYPE(bar_type) :: barra
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  CALL io_push_title('Compute nac_vec of drhox2')
  !
  band_group_myoffset = band_group%myoffset
  !
  IF(nspin == 2) THEN
     factor = 1._DP
  ELSE
     factor = 0.5_DP
  ENDIF
  !
#if defined(__CUDA)
  CALL allocate_forces_gpu()
#endif
  !
  ALLOCATE(nacvec_drhox2(n))
  ALLOCATE(nacveclc(3, nat))
  ALLOCATE(rdrhox2(dffts%nnr, nspin))
  ALLOCATE(dvpsi(npwx, band_group%nlocx, 3))
  ALLOCATE(aux1(npwx, band_group%nlocx))
  ALLOCATE(aux2(npwx, band_group%nlocx))
  !$acc enter data create(dvpsi,aux1,aux2)
  !
  nacvec_drhox2(:) = 0._DP
  !
  CALL start_bar_type(barra,'f_drhox2',kpt_pool%nloc*nat)
  !
  ! nonlocal part
  !
  DO iks = 1,kpt_pool%nloc
     !
     IF(l_spin_flip) THEN
        iks_do = flks(iks)
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
     !
     CALL g2_kin(iks)
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
#if defined(__CUDA)
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.TRUE.)
#else
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.FALSE.)
#endif
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     this_wk = wk(iks)*factor
     !
     ! ... read in GS wavefunctions iks
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     !$acc parallel loop collapse(2) present(aux1,evc)
     DO lbnd = 1,nbnd_do
        !
        ! ibnd = band_group%l2g(lbnd)+n_trunc_bands
        !
        DO ig = 1,npw
           ibnd = band_group_myoffset+lbnd+n_trunc_bands
           aux1(ig,lbnd) = evc(ig,ibnd)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     DO ia = 1,nat
        !
        ! 1) | dvpsi_i >
        !
        CALL wbse_get_dvpsi_gamma_nonlocal(ia, aux1, dvpsi)
        !
        ! 2) nacvec_drhox2 = < evc_iv2 | dvpsi_ia_iv >
        !
        !$acc host_data use_device(evc,dvgdvg_mat,aux2)
        CALL DGEMM('N', 'N', 2*npw, band_group%nloc, nbndval-n_trunc_bands, 1._DP, &
        & evc(1,n_trunc_bands+1), 2*npwx, dvgdvg_mat(1,1,iks), nbndval0x-n_trunc_bands, &
        & 0._DP, aux2, 2*npwx)
        !$acc end host_data
        !
        DO ipol = 1,3
           !
           reduce = 0._DP
           !
           !$acc parallel loop collapse(2) reduction(+:reduce) present(aux2,dvpsi) copy(reduce)
           DO lbnd = 1,nbnd_do
              DO ig = 1,npw
                 reduce = reduce + REAL(aux2(ig,lbnd),KIND=DP)*REAL(dvpsi(ig,lbnd,ipol),KIND=DP) &
                 &               + AIMAG(aux2(ig,lbnd))*AIMAG(dvpsi(ig,lbnd,ipol))
              ENDDO
           ENDDO
           !$acc end parallel
           !
           reduce = 2._DP*reduce
           !
           IF(gstart == 2) THEN
              !$acc parallel loop reduction(+:reduce) present(aux2,dvpsi) copy(reduce)
              DO lbnd = 1,nbnd_do
                 reduce = reduce - REAL(aux2(1,lbnd),KIND=DP)*REAL(dvpsi(1,lbnd,ipol),KIND=DP)
              ENDDO
              !$acc end parallel
           ENDIF
           !
           nacvec_drhox2(3*ia-3+ipol) = nacvec_drhox2(3*ia-3+ipol) - this_wk*reduce
           !
        ENDDO
        !
        CALL update_bar_type(barra,'f_drhox2',1)
        !
     ENDDO
     !
  ENDDO
  !
  CALL mp_sum(nacvec_drhox2,intra_bgrp_comm)
  CALL mp_sum(nacvec_drhox2,inter_bgrp_comm)
  CALL mp_sum(nacvec_drhox2,inter_pool_comm)
  !
  CALL stop_bar_type(barra,'f_drhox2')
  !
  ! local part
  !
  rdrhox2(:,:) = REAL(drhox2,KIND=DP)
  !
  IF(nspin == 2) THEN
     rdrhox2(:,1) = rdrhox2(:,1)+rdrhox2(:,2)
  ENDIF
  !
  CALL force_lc(nat, tau, ityp, ntyp, alat, omega, ngm, ngl, igtongl, g, rdrhox2(:,1), gstart, &
  & gamma_only, vloc, nacveclc)
  !
  nacveclc(:,:) = -factor*nacveclc
  !
  DO ia = 1,nat
     DO ipol = 1,3
        nacvec_drhox2(3*ia-3+ipol) = ( nacvec_drhox2(3*ia-3+ipol)+nacveclc(ipol,ia) ) / (-omega_JI)
     ENDDO
  ENDDO
  !
  nac_vec(:) = nac_vec+nacvec_drhox2
  !
  CALL io_push_title('eenac drhox2')
  !
  DO ia = 1,nat
     !
     WRITE(stdout, 9035) ia, ityp(ia), (nacvec_drhox2(3*ia-3+ipol), ipol = 1,3)
     !
  ENDDO
  !
  WRITE(stdout,*)
  !
  IF(mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('output.nac_vec.eenac_drhox2', nacvec_drhox2(1:n))
     !
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     !
     CALL json%destroy()
     !
  ENDIF
  !
#if defined(__CUDA)
  CALL deallocate_forces_gpu()
#endif
  !
  DEALLOCATE(nacvec_drhox2)
  DEALLOCATE(nacveclc)
  DEALLOCATE(rdrhox2)
  !$acc exit data delete(dvpsi,aux1,aux2)
  DEALLOCATE(dvpsi)
  DEALLOCATE(aux1)
  DEALLOCATE(aux2)
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   eenac = ',3F14.8)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_nacvec_drhoz_nac(n, zvector, nac_vec)
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ntyp=>nsp,ityp,tau
  USE cell_base,            ONLY : alat,omega
  USE gvect,                ONLY : g,gstart,ngm,ngl,igtongl
  USE uspp,                 ONLY : nkb,vkb
  USE uspp_init,            ONLY : init_us_2
  USE pwcom,                ONLY : isk,igk_k,lsda,current_spin,nspin,current_k,ngk,npwx,npw,xk,wk
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE westcom,              ONLY : iuwfc,lrwfc,logfile,nbnd_occ,n_trunc_bands,l_spin_flip,l_genac,&
                                 & l_eenac,do_eenac
  USE vlocal,               ONLY : vloc
  USE control_flags,        ONLY : gamma_only
  USE distribution_center,  ONLY : kpt_pool,band_group
  USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_pool_comm,inter_bgrp_comm,&
                                 & intra_bgrp_comm
  USE json_module,          ONLY : json_file
  USE mp_world,             ONLY : mpime,root
  USE io_push,              ONLY : io_push_title
  USE wavefunctions,        ONLY : evc
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_forces_gpu,deallocate_forces_gpu
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: n
  COMPLEX(DP), INTENT(IN) :: zvector(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  REAL(DP), INTENT(INOUT) :: nac_vec(n)
  !
  ! Workspace
  !
  COMPLEX(DP), ALLOCATABLE :: dvpsi(:,:,:), aux1(:,:), drhoz(:,:)
  INTEGER :: iks, iks_do, nbndval, nbnd_do, ia, ipol, lbnd, ibnd, ig
  INTEGER :: band_group_myoffset
  REAL(DP) :: reduce, factor, this_wk
  REAL(DP), ALLOCATABLE :: nacvec_drhoz(:), nacveclc(:,:), rdrhoz(:,:)
  TYPE(json_file) :: json
  INTEGER :: iunit
  TYPE(bar_type) :: barra
  CHARACTER(LEN=5) :: label
  !
  IF(l_genac .AND. .NOT. do_eenac) THEN
     label = 'genac'
  ELSEIF(l_eenac .AND. do_eenac) THEN
     label = 'eenac'
  ELSE
     CALL errore('wbse_nacvec_drhoz_nac','unexpected error',1)
  ENDIF
  !
  CALL io_push_title('Compute nac_vec of Z vector')
  !
  band_group_myoffset = band_group%myoffset
  !
  IF(nspin == 2) THEN
     factor = 1._DP
  ELSE
     factor = 0.5_DP
  ENDIF
  !
#if defined(__CUDA)
  CALL allocate_forces_gpu()
#endif
  !
  ALLOCATE(nacvec_drhoz(n))
  ALLOCATE(nacveclc(3, nat))
  ALLOCATE(rdrhoz(dffts%nnr, nspin))
  ALLOCATE(dvpsi(npwx, band_group%nlocx, 3))
  ALLOCATE(drhoz(dffts%nnr, nspin))
  ALLOCATE(aux1(npwx, band_group%nlocx))
  !$acc enter data create(dvpsi,drhoz,aux1)
  !
  nacvec_drhoz(:) = 0._DP
  !
  CALL start_bar_type(barra,'f_drhoxz',kpt_pool%nloc*nat)
  !
  ! nonlocal part
  !
  DO iks = 1,kpt_pool%nloc
     !
     ! Z vector always spin-conserving
     !
     IF(l_spin_flip) THEN
        iks_do = iks
     ELSE
        iks_do = iks
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF(lsda) current_spin = isk(iks)
     !
     CALL g2_kin(iks)
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
#if defined(__CUDA)
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.TRUE.)
#else
     IF(nkb > 0) CALL init_us_2(ngk(iks),igk_k(1,iks),xk(1,iks),vkb,.FALSE.)
#endif
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(iks)
     this_wk = wk(iks)*factor
     !
     ! ... read in GS wavefunctions iks
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     !$acc parallel loop collapse(2) present(aux1,evc)
     DO lbnd = 1,nbnd_do
        !
        ! ibnd = band_group%l2g(lbnd)+n_trunc_bands
        !
        DO ig = 1,npw
           ibnd = band_group_myoffset+lbnd+n_trunc_bands
           aux1(ig,lbnd) = evc(ig,ibnd)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     DO ia = 1,nat
        !
        ! 1) | dvpsi_i >
        !
        CALL wbse_get_dvpsi_gamma_nonlocal(ia, aux1, dvpsi)
        !
        ! 2) nacvec_drhoz_i = < z_vector | dvpsi_i >
        !
        DO ipol = 1,3
           !
           reduce = 0._DP
           !
           !$acc parallel loop collapse(2) reduction(+:reduce) present(zvector,dvpsi) copy(reduce)
           DO lbnd = 1,nbnd_do
              DO ig = 1,npw
                 reduce = reduce &
                 & + REAL(zvector(ig,lbnd,iks),KIND=DP)*REAL(dvpsi(ig,lbnd,ipol),KIND=DP) &
                 & + AIMAG(zvector(ig,lbnd,iks))*AIMAG(dvpsi(ig,lbnd,ipol))
              ENDDO
           ENDDO
           !$acc end parallel
           !
           reduce = 2._DP*reduce
           !
           IF(gstart == 2) THEN
              !$acc parallel loop reduction(+:reduce) present(zvector,dvpsi) copy(reduce)
              DO lbnd = 1,nbnd_do
                 reduce = reduce - REAL(zvector(1,lbnd,iks),KIND=DP)*REAL(dvpsi(1,lbnd,ipol),KIND=DP)
              ENDDO
             !$acc end parallel
           ENDIF
           !
           nacvec_drhoz(3*ia-3+ipol) = nacvec_drhoz(3*ia-3+ipol) + this_wk*reduce
           !
        ENDDO
        !
        CALL update_bar_type(barra,'f_drhoxz',1)
        !
     ENDDO
     !
  ENDDO
  !
  CALL mp_sum(nacvec_drhoz,intra_bgrp_comm)
  CALL mp_sum(nacvec_drhoz,inter_bgrp_comm)
  CALL mp_sum(nacvec_drhoz,inter_pool_comm)
  !
  CALL stop_bar_type(barra,'f_drhoxz')
  !
  ! local part
  !
  CALL wbse_calc_dens(zvector, drhoz, .FALSE.)
  !
  rdrhoz(:,:) = REAL(drhoz,KIND=DP)
  !
  IF(nspin == 2) THEN
     rdrhoz(:,1) = rdrhoz(:,1)+rdrhoz(:,2)
  ENDIF
  !
  CALL force_lc(nat, tau, ityp, ntyp, alat, omega, ngm, ngl, igtongl, g, rdrhoz(:,1), gstart, &
  & gamma_only, vloc, nacveclc)
  !
  nacveclc(:,:) = -factor*nacveclc
  !
  DO ia = 1,nat
     DO ipol = 1,3
        nacvec_drhoz(3*ia-3+ipol) = nacvec_drhoz(3*ia-3+ipol) + nacveclc(ipol,ia)
     ENDDO
  ENDDO
  !
  nac_vec(:) = nac_vec+nacvec_drhoz
  !
  CALL io_push_title(label//' drhoz')
  !
  DO ia = 1,nat
     !
     WRITE(stdout, 9035) ia, ityp(ia), label, (nacvec_drhoz(3*ia-3+ipol), ipol = 1,3)
     !
  ENDDO
  !
  WRITE(stdout,*)
  !
  IF(mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('output.nac_vec.'//label//'_drhoz', nacvec_drhoz(1:n))
     !
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     !
     CALL json%destroy()
     !
  ENDIF
  !
#if defined(__CUDA)
  CALL deallocate_forces_gpu()
#endif
  !
  DEALLOCATE(nacvec_drhoz)
  DEALLOCATE(nacveclc)
  DEALLOCATE(rdrhoz)
  !$acc exit data delete(dvpsi,drhoz,aux1)
  DEALLOCATE(dvpsi)
  DEALLOCATE(drhoz)
  DEALLOCATE(aux1)
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   ',A,' = ',3F14.8)
  !
END SUBROUTINE
!
END MODULE
