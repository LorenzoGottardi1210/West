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
! 
!-----------------------------------------------------------------------
SUBROUTINE wbse_calc_genac(dvg_exc_tmp)
  !-----------------------------------------------------------------------
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ityp
  USE pwcom,                ONLY : nspin,npwx
  USE noncollin_module,     ONLY : npol
  USE fft_base,             ONLY : dffts
  USE westcom,              ONLY : logfile,nbndval0x,n_trunc_bands
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
  COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol, band_group%nlocx, kpt_pool%nloc)
  !
  ! Workspace
  !
  INTEGER :: iks, n, ia, ipol
  REAL(DP), ALLOCATABLE :: nac_vec(:)
  REAL(DP) :: sumnac_vec
  COMPLEX(DP), ALLOCATABLE :: z_rhs_vec(:,:,:), zvector(:,:,:)
  TYPE(json_file) :: json
  INTEGER :: iunit
  !
  CALL start_clock('calc_geNAC')
  !
  CALL io_push_title('Compute geNAC')
  !
  n = 3 * nat
  !
  ALLOCATE(nac_vec(n))
  nac_vec(:) = 0._DP
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
  !!! SPV
  z_rhs_vec = dvg_exc_tmp
  CALL solve_zvector_eq_cg(z_rhs_vec, zvector)
  !!!
  !
#if defined(__CUDA)
  CALL deallocate_bse_gpu()
#endif
  !
  CALL wbse_nacvec_drhoz_genac(n, zvector, nac_vec)
  !
  !$acc exit data delete(z_rhs_vec,zvector)
  DEALLOCATE(z_rhs_vec)
  DEALLOCATE(zvector)
  !
  CALL io_push_title('geNAC total')
  !
  DO ia = 1,nat
     !
     WRITE(stdout, 9035) ia, ityp(ia), (nac_vec(3*ia-3+ipol), ipol = 1,3)
     !
  ENDDO
  !
  IF(mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('output.nac_vec.genac_total', nac_vec(1:n)) !!! SPV
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
  CALL io_push_title('geNAC corrected')
  !
  DO ia = 1,nat
     !
     WRITE(stdout, 9035) ia, ityp(ia), (nac_vec(3*ia-3+ipol), ipol=1,3)
     !
  ENDDO
  !
  WRITE(stdout,*)
  !
  IF(mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('output.nac_vec.genac_corrected', nac_vec(1:n))
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
  !
  CALL stop_clock('calc_geNAC')
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   geNAC = ',3F14.8)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_nacvec_drhoz_genac(n, zvector, nac_vec)
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
  USE westcom,              ONLY : iuwfc,lrwfc,logfile,nbnd_occ,n_trunc_bands,l_spin_flip
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
  !
  CALL io_push_title('Compute geNAC of Z vector')
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
        CALL wbse_get_dvpsi_gamma_nonlocal_genac(ia, aux1, dvpsi)
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
           !!! SPV no need for the c.c. in the GE-NAC
           ! nacvec_drhoz(3*ia-3+ipol) = nacvec_drhoz(3*ia-3+ipol) + 2._DP*this_wk*reduce
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
  !!! SPV no need for the c.c. in the GE-NAC
  ! drhoz(:,:) = 2._DP*drhoz
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
  CALL io_push_title('geNAC drhoz')
  !
  DO ia = 1,nat
     !
     WRITE(stdout, 9035) ia, ityp(ia), (nacvec_drhoz(3*ia-3+ipol), ipol = 1,3)
     !
  ENDDO
  !
  WRITE(stdout,*)
  !
  IF(mpime == root) THEN
     !
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
     CALL json%add('output.nac_vec.genac_drhoz', nacvec_drhoz(1:n))
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
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   geNAC = ',3F14.8)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_get_dvpsi_gamma_nonlocal_genac(i_at, dvg_tmp, dvpsi)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat,ityp,ntyp=>nsp
  USE cell_base,            ONLY : tpiba
  USE fft_interfaces,       ONLY : fwfft,invfft
  USE gvect,                ONLY : g,gstart
  USE noncollin_module,     ONLY : npol
  USE uspp_param,           ONLY : nh
  USE uspp,                 ONLY : dvan,vkb
  USE pwcom,                ONLY : npw,npwx
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE distribution_center,  ONLY : band_group
#if defined(__CUDA)
  USE cublas
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: i_at
  COMPLEX(DP), INTENT(IN) :: dvg_tmp(npwx*npol, band_group%nlocx)
  COMPLEX(DP), INTENT(OUT) :: dvpsi(npwx, band_group%nlocx, 3)
  !
  ! Workspace
  !
  INTEGER :: ia, ib, ig, nt, ih, jkb, ic, nh_nt
  INTEGER :: band_group_nloc
  COMPLEX(DP) :: factor
  REAL(DP), ALLOCATABLE :: bec1(:,:), bec2(:,:)
  COMPLEX(DP), ALLOCATABLE :: work(:,:)
  !
  !$acc kernels present(dvpsi)
  dvpsi(:,:,:) = (0._DP,0._DP)
  !$acc end kernels
  !
  band_group_nloc = band_group%nloc
  factor = tpiba*(0._DP,-1._DP)
  !
  jkb = 0
  DO nt = 1,ntyp
     nh_nt = nh(nt)
     DO ia = 1,nat
        IF(ityp(ia) == nt) THEN
           IF(ia == i_at) EXIT
           jkb = jkb+nh_nt
        ENDIF
     ENDDO
     IF(ia == i_at) EXIT
  ENDDO
  !
  IF(nh_nt < 1) RETURN
  !
  ALLOCATE(work(npwx,nh_nt))
  ALLOCATE(bec1(nh_nt,band_group%nlocx))
  ALLOCATE(bec2(nh_nt,band_group%nlocx))
  !$acc enter data create(work,bec1,bec2)
  !
  DO ic = 1,3
     !
     ! first term: sum_l sum_G' [ i V_l(G) V^*_l(G') (G'*u) psi(G')
     !
     !$acc parallel loop collapse(2) present(work,vkb,g)
     DO ih = 1,nh_nt
        DO ig = 1,npw
           work(ig,ih) = vkb(ig,jkb+ih)*g(ic,ig)*factor
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc host_data use_device(work,dvg_tmp,bec1)
     CALL DGEMM('C', 'N', nh_nt, band_group%nloc, 2*npw, 2._DP, work, 2*npwx, dvg_tmp, 2*npwx, &
     & 0._DP, bec1, nh_nt)
     !$acc end host_data
     !
     IF(gstart == 2) THEN
        !$acc parallel loop collapse(2) present(bec1,work,dvg_tmp)
        DO ib = 1,band_group_nloc
           DO ih = 1,nh_nt
              bec1(ih,ib) = bec1(ih,ib) - work(1,ih)*dvg_tmp(1,ib)
           ENDDO
        ENDDO
        !$acc end parallel
     ENDIF
     !
     !$acc host_data use_device(bec1)
     CALL mp_sum(bec1,intra_bgrp_comm)
     !$acc end host_data
     !
     !$acc parallel loop collapse(2) present(bec1,dvan)
     DO ib = 1,band_group_nloc
        DO ih = 1,nh_nt
           bec1(ih,ib) = dvan(ih,ih,nt)*bec1(ih,ib)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc host_data use_device(vkb,bec1,dvpsi)
     CALL DGEMM('N', 'N', 2*npw, band_group%nloc, nh_nt, 1._DP, vkb(1,jkb+1), 2*npwx, bec1, nh_nt, &
     & 1._DP, dvpsi(1,1,ic), 2*npwx)
     !$acc end host_data
     !
     ! second term: sum_l sum_G' [-i (G*u) V_l(G) V^*_l(G') psi(G')
     !
     !$acc host_data use_device(vkb,dvg_tmp,bec2)
     CALL DGEMM('C', 'N', nh_nt, band_group%nloc, 2*npw, 2._DP, vkb(1,jkb+1), 2*npwx, dvg_tmp, &
     & 2*npwx, 0._DP, bec2, nh_nt)
     !$acc end host_data
     !
     IF(gstart == 2) THEN
        !$acc parallel loop collapse(2) present(bec2,vkb,dvg_tmp)
        DO ib = 1,band_group_nloc
           DO ih = 1,nh_nt
              bec2(ih,ib) = bec2(ih,ib) - vkb(1,jkb+ih)*dvg_tmp(1,ib)
           ENDDO
        ENDDO
        !$acc end parallel
     ENDIF
     !
     !$acc host_data use_device(bec2)
     CALL mp_sum(bec2,intra_bgrp_comm)
     !$acc end host_data
     !
     !$acc parallel loop collapse(2) present(bec2,dvan)
     DO ib = 1,band_group_nloc
        DO ih = 1,nh_nt
           bec2(ih,ib) = dvan(ih,ih,nt)*bec2(ih,ib)
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc host_data use_device(work,bec2,dvpsi)
     CALL DGEMM('N', 'N', 2*npw, band_group%nloc, nh_nt, 1._DP, work, 2*npwx, bec2, nh_nt, 1._DP, &
     & dvpsi(1,1,ic), 2*npwx)
     !$acc end host_data
     !
  ENDDO
  !
  !$acc exit data delete(work,bec1,bec2)
  DEALLOCATE(work)
  DEALLOCATE(bec1)
  DEALLOCATE(bec2)
  !
END SUBROUTINE
