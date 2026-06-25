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
! Yu Jin, Stefano Paolo Villani, Victor Yu
!
!-----------------------------------------------------------------------
MODULE rhs_zvector
  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: build_rhs_zvector_eq
  PUBLIC :: build_rhs_zvector_eq_eenac
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE build_rhs_zvector_eq(dvg_exc_tmp,dvgdvg_mat,drhox1,drhox2,z_rhs_vec)
      !-----------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE io_push,              ONLY : io_push_title
       USE pwcom,                ONLY : ngk,npw,npwx,nspin
      USE noncollin_module,     ONLY : npol
      USE fft_base,             ONLY : dffts
      USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_bse,&
                                      & l_hybrid_tddft,l_bse_triplet
      USE mp,                   ONLY : mp_bcast
      USE buffers,              ONLY : get_buffer
      USE mp_global,            ONLY : inter_image_comm,my_image_id
      USE wavefunctions,        ONLY : evc
      USE distribution_center,  ONLY : kpt_pool,band_group
#if defined(__CUDA)
      USE west_gpu,             ONLY : reallocate_ps_gpu
#endif
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands,band_group%nlocx,kpt_pool%nloc)
      COMPLEX(DP), INTENT(IN) :: drhox1(dffts%nnr,nspin),drhox2(dffts%nnr,nspin)
      COMPLEX(DP), INTENT(OUT) :: z_rhs_vec(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      !
      ! Workspace
      !
      INTEGER :: iks,lbnd,ibnd,nbndval,nbnd_do
      !
      CALL start_clock('build_zvec')
      !
      CALL io_push_title('Build the RHS of the Z vector equation')
      !
      !$acc kernels present(z_rhs_vec)
      z_rhs_vec = (0._DP,0._DP)
      !$acc end kernels
      !
      ! part1: moved to the end
      !
      ! part2: d < a | K1e | a > / d | v >
      !
      IF(.NOT. l_bse_triplet) CALL rhs_zvector_part2(dvg_exc_tmp,z_rhs_vec,.FALSE.)
      !
      ! part3: d^2 vxc / d rho^2 contribution to d < a | K1e | a > / d | v >
      !
      IF((.NOT. l_bse) .AND. (.NOT. l_bse_triplet)) &
      & CALL rhs_zvector_part3(dvg_exc_tmp,z_rhs_vec,.FALSE.)
      !
      ! part4: d < a | K1d | a > / d | v >
      !
      IF(l_hybrid_tddft .OR. l_bse) CALL rhs_zvector_part4(dvg_exc_tmp,z_rhs_vec,.FALSE.)
      !
      ! part1: d < a | D | a > / d | v >
      !
      CALL rhs_zvector_part1(dvg_exc_tmp,dvgdvg_mat,drhox1,drhox2,z_rhs_vec,.FALSE.)
      !
      DO iks = 1,kpt_pool%nloc
         !
         nbndval = nbnd_occ(iks)
         !
         nbnd_do = 0
         DO lbnd = 1,band_group%nloc
            ibnd = band_group%l2g(lbnd)+n_trunc_bands
            IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
         ENDDO
         !
         npw = ngk(iks)
         !
         ! ... read in GS wavefunctions iks
         !
         IF(kpt_pool%nloc > 1) THEN
            IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
            CALL mp_bcast(evc,0,inter_image_comm)
            !$acc update device(evc)
         ENDIF
         !
         ! Pc[k]*z_rhs_vec
         !
#if defined(__CUDA)
         CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
         !
         CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,z_rhs_vec(:,:,iks),(1._DP,0._DP))
         !
      ENDDO
      !
      CALL stop_clock('build_zvec')
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE build_rhs_zvector_eq_eenac(dvg_exc_tmp_I,dvg_exc_tmp_J,dvgdvg_mat,dvgdvg_mat_JI,drhox1,drhox2,z_rhs_vec,omega_JI)
      !-----------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE io_push,              ONLY : io_push_title
      USE pwcom,                ONLY : ngk,npw,npwx,nspin
      USE noncollin_module,     ONLY : npol
      USE fft_base,             ONLY : dffts
      USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_bse,&
                                     & l_hybrid_tddft,l_bse_triplet
      USE mp,                   ONLY : mp_bcast
      USE buffers,              ONLY : get_buffer
      USE mp_global,            ONLY : inter_image_comm,my_image_id
      USE wavefunctions,        ONLY : evc
      USE distribution_center,  ONLY : kpt_pool,band_group
#if defined(__CUDA)
      USE west_gpu,             ONLY : reallocate_ps_gpu
#endif
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp_I(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp_J(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      REAL(DP), INTENT(IN) :: omega_JI
      REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands,band_group%nlocx,kpt_pool%nloc)
      REAL(DP), INTENT(IN) :: dvgdvg_mat_JI(nbndval0x-n_trunc_bands,band_group%nlocx,kpt_pool%nloc)
      COMPLEX(DP), INTENT(IN) :: drhox1(dffts%nnr,nspin),drhox2(dffts%nnr,nspin)
      COMPLEX(DP), INTENT(OUT) :: z_rhs_vec(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      !
      ! Workspace
      !
      INTEGER :: iks,lbnd,ibnd,nbndval,nbnd_do
      !
      CALL start_clock('build_zvec')
      !
      CALL io_push_title('Build the RHS of the Z vector equation')
      !
      !$acc kernels present(z_rhs_vec)
      z_rhs_vec = (0._DP,0._DP)
      !$acc end kernels
      !
      ! part1: moved to the end
      !
      ! part2: d < a | K1e | a > / d | v >
      !
      IF(.NOT. l_bse_triplet) THEN
         !
         ! Two calls because of derivative wrt real and complex orbitals
         !
         CALL rhs_zvector_part2(dvg_exc_tmp_I,z_rhs_vec,.TRUE.,dvg_exc_tmp_J)
         CALL rhs_zvector_part2(dvg_exc_tmp_J,z_rhs_vec,.TRUE.,dvg_exc_tmp_I)
      ENDIF
      !
      ! part3: d^2 vxc / d rho^2 contribution to d < a | K1e | a > / d | v >
      !
      IF((.NOT. l_bse) .AND. (.NOT. l_bse_triplet)) &
      & CALL rhs_zvector_part3(dvg_exc_tmp_I,z_rhs_vec,.TRUE.,dvg_exc_tmp_J)
      !
      ! part4: d < a | K1d | a > / d | v >
      !
      IF(l_hybrid_tddft) THEN
         CALL rhs_zvector_part4(dvg_exc_tmp_I,z_rhs_vec,.TRUE.,dvg_exc_tmp_J)
      ELSEIF(l_bse) THEN
         CALL errore('build_rhs_zvector_eq_eenac','BSE NACs not implemented',1)
      ENDIF
      !
      ! part1: d < a | D | a > / d | v >
      !
      CALL rhs_zvector_part1(dvg_exc_tmp_I,dvgdvg_mat,drhox1,drhox2,z_rhs_vec,.TRUE.,dvg_exc_tmp_J,&
      & dvgdvg_mat_JI)
      !
      DO iks = 1,kpt_pool%nloc
         !
         nbndval = nbnd_occ(iks)
         !
         nbnd_do = 0
         DO lbnd = 1,band_group%nloc
            ibnd = band_group%l2g(lbnd)+n_trunc_bands
            IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
         ENDDO
         !
         npw = ngk(iks)
         !
         ! ... read in GS wavefunctions iks
         !
         IF(kpt_pool%nloc > 1) THEN
            IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
            CALL mp_bcast(evc,0,inter_image_comm)
            !$acc update device(evc)
         ENDIF
         !
         ! Pc[k]*z_rhs_vec
         !
#if defined(__CUDA)
         CALL reallocate_ps_gpu(nbndval,nbnd_do)
#endif
         !
         CALL apply_alpha_pc_to_m_wfcs(nbndval,nbnd_do,z_rhs_vec(:,:,iks),(1._DP,0._DP))
         !
      ENDDO
      !
      !$acc kernels present(z_rhs_vec)
      z_rhs_vec(:,:,:) = -z_rhs_vec / omega_JI
      !$acc end kernels
      !
      CALL stop_clock('build_zvec')
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE rhs_zvector_part1(dvg_exc_tmp,dvgdvg_mat,drhox1,drhox2,z_rhs_vec,l_nac,dvg_exc_tmp_J,dvgdvg_mat_JI)
      !-----------------------------------------------------------------------
      !
      USE io_global,            ONLY : stdout
      USE kinds,                ONLY : DP
      USE io_push,              ONLY : io_push_title
      USE gvect,                ONLY : gstart
      USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_bse,&
                                     & l_hybrid_tddft,l_spin_flip,evc1_all,evc1J_all
      USE pwcom,                ONLY : isk,lsda,nspin,current_spin,current_k,ngk,npwx,npw
      USE mp,                   ONLY : mp_bcast
      USE buffers,              ONLY : get_buffer
      USE noncollin_module,     ONLY : npol
      USE fft_base,             ONLY : dffts
      USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,&
                                     & double_invfft_gamma
      USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
      USE distribution_center,  ONLY : kpt_pool,band_group
      USE mp_global,            ONLY : inter_image_comm,my_image_id
      USE wbse_dv,              ONLY : wbse_dv_setup,wbse_dv_of_drho
      USE wbse_bgrp,            ONLY : gather_bands
      USE west_mp,              ONLY : west_mp_wait
      USE xc_lib,               ONLY : xclib_dft_is
      USE wavefunctions,        ONLY : evc,psic
#if defined(__CUDA)
      USE cublas
#endif
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      REAL(DP), INTENT(IN) :: dvgdvg_mat(nbndval0x-n_trunc_bands,band_group%nlocx,kpt_pool%nloc)
      COMPLEX(DP), INTENT(IN) :: drhox1(dffts%nnr,nspin),drhox2(dffts%nnr,nspin)
      COMPLEX(DP), INTENT(INOUT) :: z_rhs_vec(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      LOGICAL, INTENT(IN) :: l_nac
      COMPLEX(DP), INTENT(IN), OPTIONAL :: dvg_exc_tmp_J(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      REAL(DP), INTENT(IN), OPTIONAL :: dvgdvg_mat_JI(nbndval0x-n_trunc_bands,band_group%nlocx,kpt_pool%nloc)
      !
      ! Workspace
      !
      INTEGER :: ibnd,jbnd,iks,iks_do,ir,ig,nbndval,nbnd_do,lbnd
      INTEGER :: dffts_nnr
      INTEGER :: req
      COMPLEX(DP), ALLOCATABLE :: dotp(:)
      COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part1(:,:,:),tmp_vec(:,:,:)
      COMPLEX(DP), ALLOCATABLE :: drhox(:,:)
      TYPE(bar_type) :: barra
      INTEGER, PARAMETER :: flks(2) = [2,1]
      !
      CALL io_push_title('Compute d <a|D|a> / d |v>')
      !
      IF(l_nac) THEN
         IF((.NOT. PRESENT(dvg_exc_tmp_J)) .OR. (.NOT. PRESENT(dvgdvg_mat_JI))) &
         & CALL errore('rhs_zvector_part1','eeNAC needs state J',1)
      ENDIF
      !
      dffts_nnr = dffts%nnr
      !
      ALLOCATE(z_rhs_vec_part1(npwx*npol,band_group%nlocx,kpt_pool%nloc))
      !$acc enter data create(z_rhs_vec_part1)
      !
      !$acc kernels present(z_rhs_vec_part1)
      z_rhs_vec_part1(:,:,:) = (0._DP,0._DP)
      !$acc end kernels
      !
      IF(xclib_dft_is('hybrid')) THEN
         !
         ALLOCATE(tmp_vec(npwx*npol,band_group%nlocx,kpt_pool%nloc))
         !$acc enter data create(tmp_vec)
         !
         !$acc kernels present(tmp_vec)
         tmp_vec(:,:,:) = (0._DP,0._DP)
         !$acc end kernels
         !
      ENDIF
      !
      ! Compute drhox
      !
      ALLOCATE(drhox(dffts%nnr,nspin))
      !
      DO iks = 1,nspin
         !
         IF(l_spin_flip) THEN
            iks_do = flks(iks)
         ELSE
            iks_do = iks
         ENDIF
         !
         drhox(:,iks) = drhox1(:,iks) + drhox2(:,iks_do)
         !
      ENDDO
      !
      !$acc enter data copyin(drhox)
      !
      CALL wbse_dv_setup(.FALSE.)
      !
      CALL wbse_dv_of_drho(drhox,.FALSE.,.FALSE.)
      !
      CALL start_bar_type(barra,'zvec1',kpt_pool%nloc)
      !
      DO iks = 1,kpt_pool%nloc
         !
         nbndval = nbnd_occ(iks)
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
         ! ... Number of G vectors for PW expansion of wfs at k
         !
         npw = ngk(iks)
         !
         ! ... read in GS wavefunctions iks
         !
         IF(kpt_pool%nloc > 1) THEN
            IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
            CALL mp_bcast(evc,0,inter_image_comm)
            !$acc update device(evc)
         ENDIF
         !
         ! ... Apply \Delta V_HXC on the z-vector
         !
         ! double bands @ gamma
         !
         DO lbnd = 1,nbnd_do-MOD(nbnd_do,2),2
            !
            ibnd = band_group%l2g(lbnd)+n_trunc_bands
            jbnd = band_group%l2g(lbnd+1)+n_trunc_bands
            !
            CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),evc(:,jbnd),psic,'Wave')
            !
            !$acc parallel loop present(psic,drhox)
            DO ir = 1,dffts_nnr
               psic(ir) = psic(ir)*CMPLX(REAL(drhox(ir,current_spin),KIND=DP),KIND=DP)
            ENDDO
            !$acc end parallel
            !
            CALL double_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part1(:,lbnd,iks),&
            & z_rhs_vec_part1(:,lbnd+1,iks),'Wave')
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
            CALL single_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),psic,'Wave')
            !
            !$acc parallel loop present(psic,drhox)
            DO ir = 1,dffts_nnr
               psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(drhox(ir,current_spin),KIND=DP),KIND=DP)
            ENDDO
            !$acc end parallel
            !
            CALL single_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part1(:,lbnd,iks),'Wave')
            !
         ENDIF
         !
         IF(l_nac) THEN
            ! Factor of 2 from the derivative wrt real and complex orbitals
            !
            !$acc kernels present(z_rhs_vec_part1)
            z_rhs_vec_part1(:,:,:) = 2._DP*z_rhs_vec_part1
            !$acc end kernels
            !
         ENDIF
         !
         IF(xclib_dft_is('hybrid')) THEN
            !
            IF(l_nac) THEN
               ! hybrid_kernel_term3 called once for a_I and once for a_J (derivative wrt real and complex orbitals)
               ! it uses global variable evc1_all and evc1J_all:
               ! first time evc1_all contains a_I and evc1J_all contains a_J, second time contents are switched
               !
               CALL gather_bands(dvg_exc_tmp(:,:,iks),evc1_all(:,:,iks),req)
               CALL west_mp_wait(req)
               CALL gather_bands(dvg_exc_tmp_J(:,:,iks),evc1J_all(:,:,iks),req)
               CALL west_mp_wait(req)
#if !defined(__GPU_MPI)
               !$acc update device(evc1_all(:,:,iks),evc1J_all(:,:,iks))
#endif
               !
               CALL hybrid_kernel_term1234(current_spin,z_rhs_vec_part1(:,:,iks),l_spin_flip,3)
               !
               ! switch the contents of evc1_all and evc1J_all
               !
               CALL gather_bands(dvg_exc_tmp_J(:,:,iks),evc1_all(:,:,iks),req)
               CALL west_mp_wait(req)
               CALL gather_bands(dvg_exc_tmp(:,:,iks),evc1J_all(:,:,iks),req)
               CALL west_mp_wait(req)
#if !defined(__GPU_MPI)
               !$acc update device(evc1_all(:,:,iks),evc1J_all(:,:,iks))
#endif
               !
               CALL hybrid_kernel_term1234(current_spin,z_rhs_vec_part1(:,:,iks),l_spin_flip,3)
               !
               ! the contents of evc1_all and evc1J_all are reverted back (may be unnecessary)
               !
               CALL gather_bands(dvg_exc_tmp(:,:,iks),evc1_all(:,:,iks),req)
               CALL west_mp_wait(req)
               CALL gather_bands(dvg_exc_tmp_J(:,:,iks),evc1J_all(:,:,iks),req)
               CALL west_mp_wait(req)
#if !defined(__GPU_MPI)
               !$acc update device(evc1_all(:,:,iks),evc1J_all(:,:,iks))
#endif
               !
            ELSE
               !
               CALL hybrid_kernel_term1234(current_spin,z_rhs_vec_part1(:,:,iks),l_spin_flip,3)
               !
            ENDIF
            !
            IF(l_spin_flip) THEN
               iks_do = flks(iks)
            ELSE
               iks_do = iks
            ENDIF
            !
            !$acc host_data use_device(evc,dvgdvg_mat,tmp_vec)
            CALL DGEMM('N','N',2*npwx*npol,nbnd_do,nbndval-n_trunc_bands,-1._DP,&
            & evc(1,1+n_trunc_bands),2*npwx*npol,dvgdvg_mat(1,1,iks_do),nbndval0x-n_trunc_bands,&
            & 0._DP,tmp_vec(1,1,iks),2*npwx*npol)
            !$acc end host_data
            !
            IF(l_nac) THEN
               !
               ! Two calls because of derivative wrt real and complex orbitals
               !
               !$acc host_data use_device(evc,dvgdvg_mat_JI,tmp_vec)
               CALL DGEMM('N','N',2*npwx*npol,nbnd_do,nbndval-n_trunc_bands,-1._DP,&
               & evc(1,1+n_trunc_bands),2*npwx*npol,dvgdvg_mat_JI(1,1,iks_do),&
               & nbndval0x-n_trunc_bands,1._DP,tmp_vec(1,1,iks),2*npwx*npol)
               !$acc end host_data
               !
            ENDIF
            !
            CALL gather_bands(tmp_vec(:,:,iks),evc1_all(:,:,iks),req)
            CALL west_mp_wait(req)
#if !defined(__GPU_MPI)
            !$acc update device(evc1_all(:,:,iks))
#endif
            !
            IF(l_hybrid_tddft) THEN
               CALL bse_kernel_gamma(current_spin,evc1_all(:,:,iks),z_rhs_vec_part1(:,:,iks),.FALSE.)
            ELSEIF(l_bse) THEN
               CALL hybrid_kernel_term1234(current_spin,z_rhs_vec_part1(:,:,iks),.FALSE.,1)
            ENDIF
            !
         ENDIF
         !
         IF(gstart == 2) THEN
            !$acc parallel loop present(z_rhs_vec_part1)
            DO lbnd = 1,nbnd_do
               z_rhs_vec_part1(1,lbnd,iks) = CMPLX(REAL(z_rhs_vec_part1(1,lbnd,iks),KIND=DP),KIND=DP)
            ENDDO
            !$acc end parallel
         ENDIF
         !
         !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part1)
         DO lbnd = 1,nbnd_do
            DO ig = 1,npw
               z_rhs_vec(ig,lbnd,iks) = z_rhs_vec(ig,lbnd,iks)-z_rhs_vec_part1(ig,lbnd,iks)
            ENDDO
         ENDDO
         !$acc end parallel
         !
         CALL update_bar_type(barra,'zvec1',1)
         !
      ENDDO
      !
      CALL stop_bar_type(barra,'zvec1')
      !
      ALLOCATE(dotp(nspin))
      !
      CALL wbse_dot(z_rhs_vec_part1,z_rhs_vec_part1,band_group%nlocx,dotp)
      !
      WRITE(stdout,*)
      WRITE(stdout,"(5x,'Norm of z_rhs_vec p1 = ',ES15.8)") SUM(REAL(dotp,KIND=DP))
      !
      DEALLOCATE(dotp)
      !$acc exit data delete(z_rhs_vec_part1)
      DEALLOCATE(z_rhs_vec_part1)
      IF(ALLOCATED(tmp_vec)) THEN
         !$acc exit data delete(tmp_vec)
         DEALLOCATE(tmp_vec)
      ENDIF
      !$acc exit data delete(drhox)
      DEALLOCATE(drhox)
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE rhs_zvector_part2(dvg_exc_tmp,z_rhs_vec,l_nac,dvg_exc_tmp_J)
      !-----------------------------------------------------------------------
      !
      USE io_global,            ONLY : stdout
      USE kinds,                ONLY : DP
      USE io_push,              ONLY : io_push_title
      USE gvect,                ONLY : gstart
      USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_bse,&
                                     & l_spin_flip,l_spin_flip_kernel
      USE pwcom,                ONLY : isk,lsda,nspin,current_spin,current_k,ngk,npwx,npw
      USE mp,                   ONLY : mp_sum,mp_bcast
      USE buffers,              ONLY : get_buffer
      USE noncollin_module,     ONLY : npol
      USE fft_base,             ONLY : dffts
      USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,&
                                     & double_invfft_gamma
      USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
      USE distribution_center,  ONLY : kpt_pool,band_group
      USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_bgrp_comm,intra_bgrp_comm
      USE wbse_dv,              ONLY : wbse_dv_of_drho,wbse_dv_of_drho_sf
      USE wavefunctions,        ONLY : evc,psic
#if defined(__CUDA)
      USE cublas
#endif
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      COMPLEX(DP), INTENT(INOUT) :: z_rhs_vec(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      LOGICAL, INTENT(IN) :: l_nac
      COMPLEX(DP), INTENT(IN), OPTIONAL :: dvg_exc_tmp_J(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      !
      ! Workspace
      !
      LOGICAL :: lrpa
      INTEGER :: ibnd,ibndp,jbnd,jbndp,kbnd,kbndp,iks,iks_do,ir,ig,nbndval,nbnd_do,lbnd,flnbndval
      INTEGER :: dffts_nnr,band_group_myoffset
      COMPLEX(DP), ALLOCATABLE :: dotp(:)
      COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part2(:,:,:),aux_g(:,:),evc_copy(:,:)
      COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
      COMPLEX(DP), ALLOCATABLE :: dpcpart(:,:)
      REAL(DP), ALLOCATABLE :: dv_vv_mat(:,:)
      TYPE(bar_type) :: barra
      INTEGER, PARAMETER :: flks(2) = [2,1]
      !
      CALL io_push_title('Compute d <a|K1e|a> / d |v>')
      !
      IF(l_nac .AND. (.NOT. PRESENT(dvg_exc_tmp_J))) &
      & CALL errore('rhs_zvector_part2','eeNAC needs state J',1)
      !
      dffts_nnr = dffts%nnr
      band_group_myoffset = band_group%myoffset
      !
      ALLOCATE(z_rhs_vec_part2(npwx*npol,band_group%nlocx,kpt_pool%nloc))
      ALLOCATE(aux_g(npwx*npol,nbndval0x-n_trunc_bands))
      ALLOCATE(dv_vv_mat(nbndval0x-n_trunc_bands,band_group%nlocx))
      ALLOCATE(dpcpart(npwx*npol,nbndval0x-n_trunc_bands))
      ALLOCATE(dvrs(dffts%nnr,nspin))
      !$acc enter data create(z_rhs_vec_part2,aux_g,dv_vv_mat,dpcpart,dvrs)
      !
      !$acc kernels present(z_rhs_vec_part2)
      z_rhs_vec_part2(:,:,:) = (0._DP,0._DP)
      !$acc end kernels
      !
      !$acc kernels present(dv_vv_mat)
      dv_vv_mat(:,:) = 0._DP
      !$acc end kernels
      !
      !$acc kernels present(dpcpart)
      dpcpart(:,:) = (0._DP,0._DP)
      !$acc end kernels
      !
      IF(.NOT. l_spin_flip) THEN
         !
         ! Calculation of the charge density response
         !
         CALL wbse_calc_dens(dvg_exc_tmp,dvrs,.FALSE.)
         !
         !$acc update device(dvrs)
         !
         lrpa = l_bse
         !
         CALL wbse_dv_of_drho(dvrs,lrpa,.FALSE.)
         !
         CALL start_bar_type(barra,'zvec2',kpt_pool%nloc)
         !
         DO iks = 1,kpt_pool%nloc
            !
            nbndval = nbnd_occ(iks)
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
            ! ... Number of G vectors for PW expansion of wfs at k
            !
            npw = ngk(iks)
            !
            ! ... read in GS wavefunctions iks
            !
            IF(kpt_pool%nloc > 1) THEN
               IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
               CALL mp_bcast(evc,0,inter_image_comm)
               !$acc update device(evc)
            ENDIF
            !
            ! ... Apply \Delta V_HXC on vector
            !
            ! double bands @ gamma
            !
            DO lbnd = 1,nbnd_do-MOD(nbnd_do,2),2
               !
               IF(l_nac) THEN
                  CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp_J(:,lbnd,iks),&
                  & dvg_exc_tmp_J(:,lbnd+1,iks),psic,'Wave')
               ELSE
                  CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),&
                  & dvg_exc_tmp(:,lbnd+1,iks),psic,'Wave')
               ENDIF
               !
               !$acc parallel loop present(psic,dvrs)
               DO ir = 1,dffts_nnr
                  psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
               ENDDO
               !$acc end parallel
               !
               CALL double_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part2(:,lbnd,iks),&
               & z_rhs_vec_part2(:,lbnd+1,iks),'Wave')
               !
            ENDDO
            !
            ! single band @ gamma
            !
            IF(MOD(nbnd_do,2) == 1) THEN
               !
               lbnd = nbnd_do
               !
               IF(l_nac) THEN
                  CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp_J(:,lbnd,iks),psic,'Wave')
               ELSE
                  CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks),psic,'Wave')
               ENDIF
               !
               !$acc parallel loop present(psic,dvrs)
               DO ir = 1,dffts_nnr
                  psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
               ENDDO
               !$acc end parallel
               !
               CALL single_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part2(:,lbnd,iks),'Wave')
               !
            ENDIF
            !
            ! Compute the second part: dv_vv_mat
            !
            ! double band @ gamma
            !
            DO jbnd = 1,(nbndval-n_trunc_bands)-MOD((nbndval-n_trunc_bands),2),2
               !
               kbnd = jbnd+1
               jbndp = jbnd+n_trunc_bands
               kbndp = kbnd+n_trunc_bands
               !
               CALL double_invfft_gamma(dffts,npw,npwx,evc(:,jbndp),evc(:,kbndp),psic,'Wave')
               !
               !$acc parallel loop present(psic,dvrs)
               DO ir = 1,dffts_nnr
                  psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
               ENDDO
               !$acc end parallel
               !
               CALL double_fwfft_gamma(dffts,npw,npwx,psic,aux_g(:,jbnd),aux_g(:,kbnd),'Wave')
               !
            ENDDO
            !
            ! single band @ gamma
            !
            IF(MOD((nbndval-n_trunc_bands),2) == 1) THEN
               !
               jbnd = nbndval-n_trunc_bands
               jbndp = jbnd+n_trunc_bands
               !
               CALL single_invfft_gamma(dffts,npw,npwx,evc(:,jbndp),psic,'Wave')
               !
               !$acc parallel loop present(psic,dvrs)
               DO ir = 1,dffts_nnr
                  psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
               ENDDO
               !$acc end parallel
               !
               CALL single_fwfft_gamma(dffts,npw,npwx,psic,aux_g(:,jbnd),'Wave')
               !
            ENDIF
            !
            ibnd = band_group%l2g(1)
            ibndp = ibnd+n_trunc_bands
            !
            CALL glbrak_gamma(aux_g,evc(:,ibndp:ibndp+nbnd_do-1),dv_vv_mat,npw,npwx,&
            & nbndval-n_trunc_bands,nbnd_do,nbndval0x-n_trunc_bands,npol)
            !
            !$acc host_data use_device(dv_vv_mat)
            CALL mp_sum(dv_vv_mat,intra_bgrp_comm)
            !$acc end host_data
            !
            IF(l_nac) THEN
               !$acc host_data use_device(dvg_exc_tmp_J,dv_vv_mat,dpcpart)
               CALL DGEMM('N','T',2*npwx*npol,nbndval-n_trunc_bands,nbnd_do,-1._DP,&
               & dvg_exc_tmp_J(1,1,iks),2*npwx*npol,dv_vv_mat,nbndval0x-n_trunc_bands,0._DP,&
               & dpcpart,2*npwx*npol)
               !$acc end host_data
            ELSE
               !$acc host_data use_device(dvg_exc_tmp,dv_vv_mat,dpcpart)
               CALL DGEMM('N','T',2*npwx*npol,nbndval-n_trunc_bands,nbnd_do,-1._DP,&
               & dvg_exc_tmp(1,1,iks),2*npwx*npol,dv_vv_mat,nbndval0x-n_trunc_bands,0._DP,dpcpart,&
               & 2*npwx*npol)
               !$acc end host_data
            ENDIF
            !
            !$acc host_data use_device(dpcpart)
            CALL mp_sum(dpcpart,inter_bgrp_comm)
            !$acc end host_data
            !
            !$acc parallel loop collapse(2) present(z_rhs_vec_part2,dpcpart)
            DO lbnd = 1,nbnd_do
               DO ig = 1,npw
                  !
                  ! ibnd = band_group%l2g(lbnd)
                  !
                  ibnd = band_group_myoffset+lbnd
                  !
                  z_rhs_vec_part2(ig,lbnd,iks) = z_rhs_vec_part2(ig,lbnd,iks)+dpcpart(ig,ibnd)
                  !
               ENDDO
            ENDDO
            !$acc end parallel
            !
            IF(gstart == 2) THEN
               !$acc parallel loop present(z_rhs_vec_part2)
               DO lbnd = 1,nbnd_do
                  z_rhs_vec_part2(1,lbnd,iks) = CMPLX(REAL(z_rhs_vec_part2(1,lbnd,iks),KIND=DP),KIND=DP)
               ENDDO
               !$acc end parallel
            ENDIF
            !
            !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part2)
            DO lbnd = 1,nbnd_do
               DO ig = 1,npw
                  z_rhs_vec(ig,lbnd,iks) = z_rhs_vec(ig,lbnd,iks)-z_rhs_vec_part2(ig,lbnd,iks)
               ENDDO
            ENDDO
            !$acc end parallel
            !
            CALL update_bar_type(barra,'zvec2',1)
            !
         ENDDO
         !
         CALL stop_bar_type(barra,'zvec2')
         !
      ELSE
         !
         IF(l_spin_flip_kernel) THEN
            !
            ALLOCATE(evc_copy(npwx,nbndval0x-n_trunc_bands))
            !$acc enter data create(evc_copy)
            !
            ! Calculation of the charge density response
            !
            CALL wbse_calc_dens(dvg_exc_tmp,dvrs,.TRUE.)
            !
            !$acc update device(dvrs)
            !
            CALL wbse_dv_of_drho_sf(dvrs)
            !
            CALL start_bar_type(barra,'zvec2',kpt_pool%nloc)
            !
            DO iks = 1,kpt_pool%nloc
               !
               iks_do = flks(iks)
               !
               nbndval = nbnd_occ(iks)
               flnbndval = nbnd_occ(iks_do)
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
               ! ... Number of G vectors for PW expansion of wfs at k
               !
               npw = ngk(iks)
               !
               ! ... read in GS wavefunctions iks
               !
               IF(kpt_pool%nloc > 1) THEN
                  IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
                  CALL mp_bcast(evc,0,inter_image_comm)
                  !$acc update device(evc)
               ENDIF
               !
               ! ... Apply \Delta V_HXC on vector
               !
               ! double bands @ gamma
               !
               DO lbnd = 1,nbnd_do-MOD(nbnd_do,2),2
                  !
                  IF(l_nac) THEN
                     CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp_J(:,lbnd,iks_do),&
                     & dvg_exc_tmp_J(:,lbnd+1,iks_do),psic,'Wave')
                  ELSE
                     CALL double_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks_do),&
                     & dvg_exc_tmp(:,lbnd+1,iks_do),psic,'Wave')
                  ENDIF
                  !
                  !$acc parallel loop present(psic,dvrs)
                  DO ir = 1,dffts_nnr
                     psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,iks_do),KIND=DP),KIND=DP)
                  ENDDO
                  !$acc end parallel
                  !
                  CALL double_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part2(:,lbnd,iks),&
                  & z_rhs_vec_part2(:,lbnd+1,iks),'Wave')
                  !
               ENDDO
               !
               ! single band @ gamma
               !
               IF(MOD(nbnd_do,2) == 1) THEN
                  !
                  lbnd = nbnd_do
                  !
                  IF(l_nac) THEN
                     CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp_J(:,lbnd,iks_do),psic,&
                     & 'Wave')
                  ELSE
                     CALL single_invfft_gamma(dffts,npw,npwx,dvg_exc_tmp(:,lbnd,iks_do),psic,'Wave')
                  ENDIF
                  !
                  !$acc parallel loop present(psic,dvrs)
                  DO ir = 1,dffts_nnr
                     psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,iks_do),KIND=DP),KIND=DP)
                  ENDDO
                  !$acc end parallel
                  !
                  CALL single_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part2(:,lbnd,iks),'Wave')
                  !
               ENDIF
               !
               !$acc parallel loop collapse(2) present(evc_copy,evc)
               DO ibnd = 1,nbndval-n_trunc_bands
                  DO ig = 1,npwx
                     evc_copy(ig,ibnd) = evc(ig,ibnd+n_trunc_bands)
                  ENDDO
               ENDDO
               !$acc end parallel
               !
               IF(kpt_pool%nloc > 1) THEN
                  IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks_do)
                  CALL mp_bcast(evc,0,inter_image_comm)
                  !$acc update device(evc)
               ENDIF
               !
               ! evc_copy -> current spin channel
               ! evc -> opposite spin channel
               !
               ! recompute nbnd_do for the opposite spin channel
               !
               nbnd_do = 0
               DO lbnd = 1,band_group%nloc
                  ibnd = band_group%l2g(lbnd)+n_trunc_bands
                  IF(ibnd > n_trunc_bands .AND. ibnd <= flnbndval) nbnd_do = nbnd_do+1
               ENDDO
               !
               ! Compute the second part: dv_vv_mat
               !
               ! double band @ gamma
               !
               DO jbnd = 1,(nbndval-n_trunc_bands)-MOD((nbndval-n_trunc_bands),2),2
                  !
                  kbnd = jbnd+1
                  !
                  CALL double_invfft_gamma(dffts,npw,npwx,evc_copy(:,jbnd),evc_copy(:,kbnd),psic,&
                  & 'Wave')
                  !
                  !$acc parallel loop present(psic,dvrs)
                  DO ir = 1,dffts_nnr
                     psic(ir) = psic(ir)*CMPLX(REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
                  ENDDO
                  !$acc end parallel
                  !
                  CALL double_fwfft_gamma(dffts,npw,npwx,psic,aux_g(:,jbnd),aux_g(:,kbnd),'Wave')
                  !
               ENDDO
               !
               ! single band @ gamma
               !
               IF(MOD((nbndval-n_trunc_bands),2) == 1) THEN
                  !
                  jbnd = nbndval-n_trunc_bands
                  !
                  CALL single_invfft_gamma(dffts,npw,npwx,evc_copy(:,jbnd),psic,'Wave')
                  !
                  !$acc parallel loop present(psic,dvrs)
                  DO ir = 1,dffts_nnr
                     psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(dvrs(ir,current_spin),KIND=DP),KIND=DP)
                  ENDDO
                  !$acc end parallel
                  !
                  CALL single_fwfft_gamma(dffts,npw,npwx,psic,aux_g(:,jbnd),'Wave')
                  !
               ENDIF
               !
               ibnd = band_group%l2g(1)
               ibndp = ibnd+n_trunc_bands
               !
               CALL glbrak_gamma(aux_g,evc(:,ibndp:ibndp+nbnd_do-1),dv_vv_mat,npw,npwx,&
               & nbndval-n_trunc_bands,nbnd_do,nbndval0x-n_trunc_bands,npol)
               !
               !$acc host_data use_device(dv_vv_mat)
               CALL mp_sum(dv_vv_mat,intra_bgrp_comm)
               !$acc end host_data
               !
               IF(l_nac) THEN
                  !$acc host_data use_device(dvg_exc_tmp_J,dv_vv_mat,dpcpart)
                  CALL DGEMM('N','T',2*npwx*npol,nbndval-n_trunc_bands,nbnd_do,-1._DP,&
                  & dvg_exc_tmp_J(1,1,iks),2*npwx*npol,dv_vv_mat,nbndval0x-n_trunc_bands,0._DP,&
                  & dpcpart,2*npwx*npol)
                  !$acc end host_data
               ELSE
                  !$acc host_data use_device(dvg_exc_tmp,dv_vv_mat,dpcpart)
                  CALL DGEMM('N','T',2*npwx*npol,nbndval-n_trunc_bands,nbnd_do,-1._DP,&
                  & dvg_exc_tmp(1,1,iks),2*npwx*npol,dv_vv_mat,nbndval0x-n_trunc_bands,0._DP,&
                  & dpcpart,2*npwx*npol)
                  !$acc end host_data
               ENDIF
               !
               !$acc host_data use_device(dpcpart)
               CALL mp_sum(dpcpart,inter_bgrp_comm)
               !$acc end host_data
               !
               ! recompute nbnd_do for the current spin channel
               !
               nbnd_do = 0
               DO lbnd = 1,band_group%nloc
                  ibnd = band_group%l2g(lbnd)+n_trunc_bands
                  IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
               ENDDO
               !
               !$acc parallel loop collapse(2) present(z_rhs_vec_part2,dpcpart)
               DO lbnd = 1,nbnd_do
                  DO ig = 1,npw
                     !
                     ! ibnd = band_group%l2g(lbnd)
                     !
                     ibnd = band_group_myoffset+lbnd
                     !
                     z_rhs_vec_part2(ig,lbnd,iks) = z_rhs_vec_part2(ig,lbnd,iks)+dpcpart(ig,ibnd)
                     !
                  ENDDO
               ENDDO
               !$acc end parallel
               !
               IF(gstart == 2) THEN
                  !$acc parallel loop present(z_rhs_vec_part2)
                  DO lbnd = 1,nbnd_do
                     z_rhs_vec_part2(1,lbnd,iks) = CMPLX(REAL(z_rhs_vec_part2(1,lbnd,iks),KIND=DP),KIND=DP)
                  ENDDO
                  !$acc end parallel
               ENDIF
               !
               !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part2)
               DO lbnd = 1,nbnd_do
                  DO ig = 1,npw
                     z_rhs_vec(ig,lbnd,iks) = z_rhs_vec(ig,lbnd,iks)-z_rhs_vec_part2(ig,lbnd,iks)
                  ENDDO
               ENDDO
               !$acc end parallel
               !
               CALL update_bar_type(barra,'zvec2',1)
               !
            ENDDO
            !
            CALL stop_bar_type(barra,'zvec2')
            !
            !$acc exit data delete(evc_copy)
            DEALLOCATE(evc_copy)
            !
         ENDIF
         !
      ENDIF
      !
      ALLOCATE(dotp(nspin))
      !
      CALL wbse_dot(z_rhs_vec_part2,z_rhs_vec_part2,band_group%nlocx,dotp)
      !
      WRITE(stdout,*)
      WRITE(stdout,"(5x,'Norm of z_rhs_vec p2 = ',ES15.8)") SUM(REAL(dotp,KIND=DP))
      !
      DEALLOCATE(dotp)
      !$acc exit data delete(z_rhs_vec_part2,aux_g,dv_vv_mat,dpcpart,dvrs)
      DEALLOCATE(z_rhs_vec_part2)
      DEALLOCATE(aux_g)
      DEALLOCATE(dv_vv_mat)
      DEALLOCATE(dpcpart)
      DEALLOCATE(dvrs)
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE rhs_zvector_part3(dvg_exc_tmp,z_rhs_vec,l_nac,dvg_exc_tmp_J)
      !-----------------------------------------------------------------------
      !
      USE io_global,            ONLY : stdout
      USE kinds,                ONLY : DP
      USE io_push,              ONLY : io_push_title
      USE gvect,                ONLY : gstart
      USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,n_trunc_bands,l_spin_flip
      USE pwcom,                ONLY : isk,lsda,nspin,current_spin,current_k,ngk,npwx,npw
      USE mp,                   ONLY : mp_bcast
      USE buffers,              ONLY : get_buffer
      USE noncollin_module,     ONLY : npol
      USE fft_base,             ONLY : dffts
      USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_fwfft_gamma,&
                                     & double_invfft_gamma
      USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
      USE distribution_center,  ONLY : kpt_pool,band_group
      USE mp_global,            ONLY : inter_image_comm,my_image_id
      USE wavefunctions,        ONLY : evc,psic
#if defined(__CUDA)
      USE cublas
#endif
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      COMPLEX(DP), INTENT(INOUT) :: z_rhs_vec(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      LOGICAL, INTENT(IN) :: l_nac
      COMPLEX(DP), INTENT(IN), OPTIONAL :: dvg_exc_tmp_J(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      !
      ! Workspace
      !
      INTEGER :: ibnd,jbnd,iks,ir,ig,nbndval,nbnd_do,lbnd
      INTEGER :: dffts_nnr
      COMPLEX(DP), ALLOCATABLE :: dotp(:)
      COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part3(:,:,:)
      COMPLEX(DP), ALLOCATABLE :: ddvxc(:,:)
      TYPE(bar_type) :: barra
      !
      CALL io_push_title('Compute d^2 vxc / d rho^2')
      !
      IF(l_nac .AND. (.NOT. PRESENT(dvg_exc_tmp_J))) &
      & CALL errore('rhs_zvector_part3','eeNAC needs state J',1)
      !
      dffts_nnr = dffts%nnr
      !
      ALLOCATE(z_rhs_vec_part3(npwx*npol,band_group%nlocx,kpt_pool%nloc))
      !$acc enter data create(z_rhs_vec_part3)
      !
      !$acc kernels present(z_rhs_vec_part3)
      z_rhs_vec_part3(:,:,:) = (0._DP,0._DP)
      !$acc end kernels
      !
      ALLOCATE(ddvxc(dffts%nnr,nspin))
      !
      IF(l_nac) THEN
         IF(.NOT. l_spin_flip) THEN
            CALL compute_ddvxc_5p_eenac(dvg_exc_tmp,dvg_exc_tmp_J,ddvxc)
         ELSE
            CALL compute_ddvxc_sf_eenac(dvg_exc_tmp,dvg_exc_tmp_J,ddvxc)
         ENDIF
      ELSE
         IF(.NOT. l_spin_flip) THEN
            CALL compute_ddvxc_5p(dvg_exc_tmp,ddvxc)
         ELSE
            CALL compute_ddvxc_sf(dvg_exc_tmp,ddvxc)
         ENDIF
      ENDIF
      !
      !$acc enter data copyin(ddvxc)
      !
      CALL start_bar_type(barra,'zvec3',kpt_pool%nloc)
      !
      DO iks = 1,kpt_pool%nloc
         !
         nbndval = nbnd_occ(iks)
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
         ! ... Number of G vectors for PW expansion of wfs at k
         !
         npw = ngk(iks)
         !
         ! ... read in GS wavefunctions iks
         !
         IF(kpt_pool%nloc > 1) THEN
            IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
            CALL mp_bcast(evc,0,inter_image_comm)
            !$acc update device(evc)
         ENDIF
         !
         ! ... Apply \Delta V_HXC
         !
         ! double bands @ gamma
         !
         DO lbnd = 1,nbnd_do-MOD(nbnd_do,2),2
            !
            ibnd = band_group%l2g(lbnd)+n_trunc_bands
            jbnd = band_group%l2g(lbnd+1)+n_trunc_bands
            !
            CALL double_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),evc(:,jbnd),psic,'Wave')
            !
            !$acc parallel loop present(psic,ddvxc)
            DO ir = 1,dffts_nnr
               psic(ir) = psic(ir)*CMPLX(REAL(ddvxc(ir,current_spin),KIND=DP),KIND=DP)
            ENDDO
            !$acc end parallel
            !
            CALL double_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part3(:,lbnd,iks),&
            & z_rhs_vec_part3(:,lbnd+1,iks),'Wave')
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
            CALL single_invfft_gamma(dffts,npw,npwx,evc(:,ibnd),psic,'Wave')
            !
            !$acc parallel loop present(psic,ddvxc)
            DO ir = 1,dffts_nnr
               psic(ir) = CMPLX(REAL(psic(ir),KIND=DP)*REAL(ddvxc(ir,current_spin),KIND=DP),KIND=DP)
            ENDDO
            !$acc end parallel
            !
            CALL single_fwfft_gamma(dffts,npw,npwx,psic,z_rhs_vec_part3(:,lbnd,iks),'Wave')
            !
         ENDIF
         !
         IF(gstart == 2) THEN
            !$acc parallel loop present(z_rhs_vec_part3)
            DO lbnd = 1,nbnd_do
               z_rhs_vec_part3(1,lbnd,iks) = CMPLX(REAL(z_rhs_vec_part3(1,lbnd,iks),KIND=DP),KIND=DP)
            ENDDO
            !$acc end parallel
         ENDIF
         !
         IF(l_nac) THEN
            !
            ! Factor of 2 from the derivative wrt real and complex orbitals
            !
            !$acc kernels present(z_rhs_vec_part3)
            z_rhs_vec_part3(:,:,:) = 2._DP*z_rhs_vec_part3
            !$acc end kernels
            !
         ENDIF
         !
         !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part3)
         DO lbnd = 1,nbnd_do
            DO ig = 1,npw
               z_rhs_vec(ig,lbnd,iks) = z_rhs_vec(ig,lbnd,iks)-z_rhs_vec_part3(ig,lbnd,iks)
            ENDDO
         ENDDO
         !$acc end parallel
         !
         CALL update_bar_type(barra,'zvec3',1)
         !
      ENDDO
      !
      CALL stop_bar_type(barra,'zvec3')
      !
      ALLOCATE(dotp(nspin))
      !
      CALL wbse_dot(z_rhs_vec_part3,z_rhs_vec_part3,band_group%nlocx,dotp)
      !
      WRITE(stdout,*)
      WRITE(stdout,"(5x,'Norm of z_rhs_vec p3 = ',ES15.8)") SUM(REAL(dotp,KIND=DP))
      !
      DEALLOCATE(dotp)
      !$acc exit data delete(z_rhs_vec_part3,ddvxc)
      DEALLOCATE(z_rhs_vec_part3)
      DEALLOCATE(ddvxc)
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_ddvxc_5p(dvg_exc_tmp,ddvxc)
      !-----------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE lsda_mod,             ONLY : nspin
      USE wvfct,                ONLY : npwx
      USE noncollin_module,     ONLY : npol
      USE fft_base,             ONLY : dffts
      USE scf,                  ONLY : rho,rho_core,rhog_core,scf_type,create_scf_type,&
                                     & destroy_scf_type
      USE distribution_center,  ONLY : kpt_pool,band_group
      USE fft_interfaces,       ONLY : fwfft
      USE westcom,              ONLY : ddvxc_fd_coeff
      USE wavefunctions,        ONLY : psic
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      COMPLEX(DP), INTENT(OUT) :: ddvxc(dffts%nnr,nspin)
      !
      ! Workspace
      !
      INTEGER :: iks,ir,indk
      REAL(DP), ALLOCATABLE :: aux_vxc(:,:,:),vxc(:,:),rdvrs(:,:)
      REAL(DP) :: etxc,vtxc
      TYPE(scf_type) :: a_rho
      COMPLEX(DP), ALLOCATABLE :: dvrs(:,:)
      !
      CALL start_clock('ddvxc_5p')
      !
      CALL create_scf_type(a_rho)
      !
      ALLOCATE(aux_vxc(dffts%nnr,nspin,5))
      ALLOCATE(vxc(dffts%nnr,nspin))
      ALLOCATE(dvrs(dffts%nnr,nspin))
      !$acc enter data create(dvrs)
      ALLOCATE(rdvrs(dffts%nnr,nspin))
      !
      ! Calculation of the charge density response
      !
      CALL wbse_calc_dens(dvg_exc_tmp,dvrs,.FALSE.)
      !
      IF(nspin == 1) THEN
         !
         rdvrs(:,1) = REAL(dvrs(:,1),KIND=DP)
         !
      ELSEIF(nspin == 2) THEN
         !
         rdvrs(:,1) = REAL(dvrs(:,1),KIND=DP) + REAL(dvrs(:,2),KIND=DP)
         rdvrs(:,2) = REAL(dvrs(:,1),KIND=DP) - REAL(dvrs(:,2),KIND=DP)
         !
      ELSEIF(nspin == 4) THEN
         !
         CALL errore('compute_ddvxc_5p','nspin == 4 not supported',1)
         !
      ENDIF
      !
      DO indk = 1,5
         !
         vxc(:,:) = 0._DP
         !
         a_rho%of_r(:,:) = rho%of_r + REAL((indk-3),KIND=DP) * ddvxc_fd_coeff * rdvrs
         !
         DO iks = 1,nspin
            !
            psic(:) = a_rho%of_r(:,iks)
            CALL fwfft('Rho',psic,dffts)
            a_rho%of_g(:,iks) = psic(dffts%nl)
            !
         ENDDO
         !
         CALL v_xc(a_rho,rho_core,rhog_core,etxc,vtxc,vxc)
         !
         aux_vxc(:,:,indk) = vxc
         !
      ENDDO
      !
      ! compute ddvxc
      !
      DO iks = 1,nspin
         DO ir = 1,dffts%nnr
            ddvxc(ir,iks) = CMPLX((-aux_vxc(ir,iks,1)+16._DP*aux_vxc(ir,iks,2) &
            &                      -30._DP*aux_vxc(ir,iks,3)+16._DP*aux_vxc(ir,iks,4) &
            &                      -aux_vxc(ir,iks,5)),KIND=DP) &
            &             / (12._DP*ddvxc_fd_coeff**2)
         ENDDO
      ENDDO
      !
      DEALLOCATE(aux_vxc)
      DEALLOCATE(vxc)
      !$acc exit data delete(dvrs)
      DEALLOCATE(dvrs)
      DEALLOCATE(rdvrs)
      !
      CALL destroy_scf_type(a_rho)
      !
      CALL stop_clock('ddvxc_5p')
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_ddvxc_sf(dvg_exc_tmp,ddvxc)
      !-----------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE lsda_mod,             ONLY : nspin
      USE wvfct,                ONLY : npwx
      USE noncollin_module,     ONLY : npol,nspin_gga
      USE xc_lib,               ONLY : xclib_dft_is
      USE fft_base,             ONLY : dffts
      USE gvect,                ONLY : g
      USE scf,                  ONLY : rho
      USE uspp,                 ONLY : nlcc_any
      USE distribution_center,  ONLY : kpt_pool,band_group
      USE westcom,              ONLY : sf_kernel,l_spin_flip_kernel,l_spin_flip_alda0,spin_flip_cut
      USE qpoint,               ONLY : xq
      USE gc_lr,                ONLY : grho,dvxc_rr,dvxc_sr,dvxc_ss,dvxc_s
      USE eqv,                  ONLY : dmuxc
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      COMPLEX(DP), INTENT(OUT) :: ddvxc(dffts%nnr,nspin)
      !
      ! Workspace
      !
      INTEGER :: ir,is,is1
      REAL(DP) :: tmp1,tmp2
      COMPLEX(DP), ALLOCATABLE :: drho_sf(:,:),drho_sf_copy(:,:)
      !
      CALL start_clock('ddvxc_sf')
      !
      IF(nlcc_any) CALL errore('compute_ddvxc_sf','nlcc_any not supported',1)
      !
      ALLOCATE(drho_sf(dffts%nnr,2))
      !$acc enter data create(drho_sf)
      ALLOCATE(drho_sf_copy(dffts%nnr,2))
      !
      CALL wbse_calc_dens(dvg_exc_tmp,drho_sf,.TRUE.)
      !
      DO ir = 1,dffts%nnr
         tmp1 = REAL(drho_sf(ir,1),KIND=DP)**2
         tmp2 = REAL(drho_sf(ir,2),KIND=DP)**2
         drho_sf_copy(ir,1) = CMPLX(tmp1+tmp2,KIND=DP)
         drho_sf_copy(ir,2) = -CMPLX(tmp1+tmp2,KIND=DP)
      ENDDO
      !
      ! divide rho_diff
      !
      DO ir = 1,dffts%nnr
         IF(ABS(rho%of_r(ir,2)) < spin_flip_cut) THEN
            drho_sf_copy(ir,1) = (0._DP,0._DP)
            drho_sf_copy(ir,2) = (0._DP,0._DP)
         ELSE
            drho_sf_copy(ir,1) = drho_sf_copy(ir,1) / rho%of_r(ir,2)
            drho_sf_copy(ir,2) = drho_sf_copy(ir,2) / rho%of_r(ir,2)
         ENDIF
      ENDDO
      !
      ddvxc(:,:) = (0._DP,0._DP)
      !
      IF(l_spin_flip_kernel) THEN
         !
         ! part 2
         !
         DO is = 1,nspin
            DO is1 = 1,nspin
               ddvxc(:,is) = ddvxc(:,is) + dmuxc(:,is,is1) * drho_sf_copy(:,is1)
            ENDDO
         ENDDO
         !
         IF(.NOT. l_spin_flip_alda0) THEN
            IF(xclib_dft_is('gradient')) THEN
               CALL dgradcorr(dffts,rho%of_r,grho,dvxc_rr,dvxc_sr,dvxc_ss,dvxc_s,xq,drho_sf_copy,&
               & nspin,nspin_gga,g,ddvxc)
            ENDIF
         ENDIF
         !
         ! part 1
         !
         !$acc update host(sf_kernel)
         !
         DO ir = 1,dffts%nnr
            ddvxc(ir,1) = ddvxc(ir,1) - sf_kernel(ir) * drho_sf_copy(ir,1)
            ddvxc(ir,2) = ddvxc(ir,2) - sf_kernel(ir) * drho_sf_copy(ir,2)
         ENDDO
         !
      ENDIF
      !
      !$acc exit data delete(drho_sf)
      DEALLOCATE(drho_sf)
      DEALLOCATE(drho_sf_copy)
      !
      CALL stop_clock('ddvxc_sf')
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_ddvxc_5p_eenac(dvg_exc_tmp,dvg_exc_tmp_J,ddvxc)
      !-----------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE lsda_mod,             ONLY : nspin
      USE wvfct,                ONLY : npwx
      USE noncollin_module,     ONLY : npol
      USE fft_base,             ONLY : dffts
      USE scf,                  ONLY : rho,rho_core,rhog_core,scf_type,create_scf_type,&
                                     & destroy_scf_type
      USE distribution_center,  ONLY : kpt_pool,band_group
      USE fft_interfaces,       ONLY : fwfft
      USE westcom,              ONLY : ddvxc_fd_coeff
      USE wavefunctions,        ONLY : psic
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp_J(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      COMPLEX(DP), INTENT(OUT) :: ddvxc(dffts%nnr,nspin)
      !
      ! Workspace
      !
      INTEGER :: iks,ir,indk
      REAL(DP), ALLOCATABLE :: aux_vxc(:,:,:),vxc(:,:),rdvrs(:,:),rdvrs_I(:,:),rdvrs_J(:,:)
      REAL(DP) :: etxc,vtxc
      TYPE(scf_type) :: a_rho
      COMPLEX(DP), ALLOCATABLE :: dvrs_I(:,:)
      COMPLEX(DP), ALLOCATABLE :: dvrs_J(:,:)
      INTEGER :: isgn
      INTEGER, PARAMETER :: signs(2) = [-1,1]
      COMPLEX(DP), ALLOCATABLE :: ddvxc_tmp(:,:,:)
      !
      CALL start_clock('ddvxc_5p')
      !
      CALL create_scf_type(a_rho)
      !
      ALLOCATE(aux_vxc(dffts%nnr,nspin,5))
      ALLOCATE(vxc(dffts%nnr,nspin))
      ALLOCATE(dvrs_I(dffts%nnr,nspin))
      ALLOCATE(dvrs_J(dffts%nnr,nspin))
      ALLOCATE(ddvxc_tmp(dffts%nnr,nspin,2))
      !$acc enter data create(dvrs_I,dvrs_J)
      ALLOCATE(rdvrs(dffts%nnr,nspin))
      ALLOCATE(rdvrs_I(dffts%nnr,nspin))
      ALLOCATE(rdvrs_J(dffts%nnr,nspin))
      !
      ! Calculation of the charge density response
      !
      CALL wbse_calc_dens(dvg_exc_tmp,dvrs_I,.FALSE.)
      CALL wbse_calc_dens(dvg_exc_tmp_J,dvrs_J,.FALSE.)
      !
      IF(nspin == 1) THEN
         !
         rdvrs_I(:,1) = REAL(dvrs_I(:,1),KIND=DP)
         rdvrs_J(:,1) = REAL(dvrs_J(:,1),KIND=DP)
         !
      ELSEIF(nspin == 2) THEN
         !
         rdvrs_I(:,1) = REAL(dvrs_I(:,1),KIND=DP) + REAL(dvrs_I(:,2),KIND=DP)
         rdvrs_I(:,2) = REAL(dvrs_I(:,1),KIND=DP) - REAL(dvrs_I(:,2),KIND=DP)
         !
         rdvrs_J(:,1) = REAL(dvrs_J(:,1),KIND=DP) + REAL(dvrs_J(:,2),KIND=DP)
         rdvrs_J(:,2) = REAL(dvrs_J(:,1),KIND=DP) - REAL(dvrs_J(:,2),KIND=DP)
         !
      ELSEIF(nspin == 4) THEN
         !
         CALL errore('compute_ddvxc_5p_eenac','nspin == 4 not supported',1)
         !
      ENDIF
      !
      DO isgn = 1,2 ! signs=[-1,+1]
         DO indk = 1,5
            !
            vxc(:,:) = 0._DP
            !
            rdvrs = rdvrs_I + signs(isgn)*rdvrs_J
            a_rho%of_r(:,:) = rho%of_r + REAL((indk-3),KIND=DP) * ddvxc_fd_coeff * rdvrs
            !
            DO iks = 1,nspin
               !
               psic(:) = a_rho%of_r(:,iks)
               CALL fwfft ('Rho',psic,dffts)
               a_rho%of_g(:,iks) = psic(dffts%nl)
               !
            ENDDO
            !
            CALL v_xc(a_rho,rho_core,rhog_core,etxc,vtxc,vxc)
            !
            aux_vxc(:,:,indk) = vxc
            !
         ENDDO
         !
         ! compute ddvxc
         !
         DO iks = 1,nspin
            DO ir = 1,dffts%nnr
               ddvxc_tmp(ir,iks,isgn) = CMPLX((-aux_vxc(ir,iks,1)+16._DP*aux_vxc(ir,iks,2) &
               &                               -30._DP*aux_vxc(ir,iks,3)+16._DP*aux_vxc(ir,iks,4) &
               &                               -aux_vxc(ir,iks,5)),KIND=DP) &
                                      / (12._DP*ddvxc_fd_coeff**2)
            ENDDO
         ENDDO
      ENDDO
      !
      ddvxc(:,:) = 0.25_DP * (ddvxc_tmp(:,:,2) - ddvxc_tmp(:,:,1))
      !
      DEALLOCATE(aux_vxc)
      DEALLOCATE(vxc)
      !$acc exit data delete(dvrs_I,dvrs_J)
      DEALLOCATE(dvrs_I)
      DEALLOCATE(dvrs_J)
      DEALLOCATE(rdvrs)
      DEALLOCATE(rdvrs_I)
      DEALLOCATE(rdvrs_J)
      !
      CALL destroy_scf_type(a_rho)
      !
      CALL stop_clock('ddvxc_5p')
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_ddvxc_sf_eenac(dvg_exc_tmp,dvg_exc_tmp_J,ddvxc)
      !-----------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE lsda_mod,             ONLY : nspin
      USE wvfct,                ONLY : npwx
      USE noncollin_module,     ONLY : npol,nspin_gga
      USE xc_lib,               ONLY : xclib_dft_is
      USE fft_base,             ONLY : dffts
      USE gvect,                ONLY : g
      USE scf,                  ONLY : rho
      USE uspp,                 ONLY : nlcc_any
      USE distribution_center,  ONLY : kpt_pool,band_group
      USE westcom,              ONLY : sf_kernel,l_spin_flip_kernel,l_spin_flip_alda0,spin_flip_cut
      USE qpoint,               ONLY : xq
      USE gc_lr,                ONLY : grho,dvxc_rr,dvxc_sr,dvxc_ss,dvxc_s
      USE eqv,                  ONLY : dmuxc
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp_J(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      COMPLEX(DP), INTENT(OUT) :: ddvxc(dffts%nnr,nspin)
      !
      ! Workspace
      !
      INTEGER :: ir,is,is1
      REAL(DP) :: tmp1,tmp2
      COMPLEX(DP), ALLOCATABLE :: drho_sf_I(:,:),drho_sf_J(:,:),drho_sf_copy(:,:)
      !
      CALL start_clock('ddvxc_sf')
      !
      IF(nlcc_any) CALL errore('compute_ddvxc_sf_eenac','nlcc_any not supported',1)
      !
      ALLOCATE(drho_sf_I(dffts%nnr,2))
      ALLOCATE(drho_sf_J(dffts%nnr,2))
      !$acc enter data create(drho_sf_I,drho_sf_J)
      ALLOCATE(drho_sf_copy(dffts%nnr,2))
      !
      CALL wbse_calc_dens(dvg_exc_tmp,drho_sf_I,.TRUE.)
      CALL wbse_calc_dens(dvg_exc_tmp_J,drho_sf_J,.TRUE.)
      !
      DO ir = 1,dffts%nnr
         tmp1 = REAL(drho_sf_I(ir,1),KIND=DP)*REAL(drho_sf_J(ir,1),KIND=DP)
         tmp2 = REAL(drho_sf_I(ir,2),KIND=DP)*REAL(drho_sf_J(ir,2),KIND=DP)
         drho_sf_copy(ir,1) = CMPLX(tmp1+tmp2,KIND=DP)
         drho_sf_copy(ir,2) = -CMPLX(tmp1+tmp2,KIND=DP)
      ENDDO
      !
      ! divide rho_diff
      !
      DO ir = 1,dffts%nnr
         IF(ABS(rho%of_r(ir,2)) < spin_flip_cut) THEN
            drho_sf_copy(ir,1) = (0._DP,0._DP)
            drho_sf_copy(ir,2) = (0._DP,0._DP)
         ELSE
            drho_sf_copy(ir,1) = drho_sf_copy(ir,1) / rho%of_r(ir,2)
            drho_sf_copy(ir,2) = drho_sf_copy(ir,2) / rho%of_r(ir,2)
         ENDIF
      ENDDO
      !
      ddvxc(:,:) = (0._DP,0._DP)
      !
      IF(l_spin_flip_kernel) THEN
         !
         ! part 2
         !
         DO is = 1,nspin
            DO is1 = 1,nspin
               ddvxc(:,is) = ddvxc(:,is) + dmuxc(:,is,is1) * drho_sf_copy(:,is1)
            ENDDO
         ENDDO
         !
         IF(.NOT. l_spin_flip_alda0) THEN
            IF(xclib_dft_is('gradient')) THEN
               CALL dgradcorr(dffts,rho%of_r,grho,dvxc_rr,dvxc_sr,dvxc_ss,dvxc_s,xq,drho_sf_copy,&
               & nspin,nspin_gga,g,ddvxc)
            ENDIF
         ENDIF
         !
         ! part 1
         !
         !$acc update host(sf_kernel)
         !
         DO ir = 1,dffts%nnr
            ddvxc(ir,1) = ddvxc(ir,1) - sf_kernel(ir) * drho_sf_copy(ir,1)
            ddvxc(ir,2) = ddvxc(ir,2) - sf_kernel(ir) * drho_sf_copy(ir,2)
         ENDDO
         !
      ENDIF
      !
      !$acc exit data delete(drho_sf_I,drho_sf_J)
      DEALLOCATE(drho_sf_I)
      DEALLOCATE(drho_sf_J)
      DEALLOCATE(drho_sf_copy)
      !
      CALL stop_clock('ddvxc_sf')
      !
    END SUBROUTINE
    !
    !-----------------------------------------------------------------------
    SUBROUTINE rhs_zvector_part4(dvg_exc_tmp,z_rhs_vec,l_nac,dvg_exc_tmp_J)
      !-----------------------------------------------------------------------
      !
      USE io_global,            ONLY : stdout
      USE kinds,                ONLY : DP
      USE io_push,              ONLY : io_push_title
      USE gvect,                ONLY : gstart
      USE westcom,              ONLY : iuwfc,lrwfc,nbnd_occ,nbndval0x,n_trunc_bands,l_bse,&
                                     & l_hybrid_tddft,l_spin_flip,evc1_all,evc1J_all
      USE pwcom,                ONLY : isk,lsda,nspin,current_spin,current_k,ngk,npwx,npw
      USE mp,                   ONLY : mp_sum,mp_bcast
      USE buffers,              ONLY : get_buffer
      USE noncollin_module,     ONLY : npol
      USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
      USE distribution_center,  ONLY : kpt_pool,band_group
      USE mp_global,            ONLY : inter_image_comm,my_image_id,inter_bgrp_comm,intra_bgrp_comm
      USE wavefunctions,        ONLY : evc
      USE wbse_bgrp,            ONLY : gather_bands
      USE west_mp,              ONLY : west_mp_wait
#if defined(__CUDA)
      USE cublas
#endif
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP), INTENT(IN) :: dvg_exc_tmp(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      COMPLEX(DP), INTENT(INOUT) :: z_rhs_vec(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      LOGICAL, INTENT(IN) :: l_nac
      COMPLEX(DP), INTENT(IN), OPTIONAL :: dvg_exc_tmp_J(npwx*npol,band_group%nlocx,kpt_pool%nloc)
      !
      ! Workspace
      !
      INTEGER :: ig,lbnd,ibnd,iks,iks_do,nbnd_do,nbndval,flnbndval
      INTEGER :: band_group_myoffset
      INTEGER :: req
      COMPLEX(DP), ALLOCATABLE :: dotp(:)
      COMPLEX(DP), ALLOCATABLE :: z_rhs_vec_part4(:,:,:),tmp_vec(:,:)
      REAL(DP), ALLOCATABLE :: dv_vv_mat(:,:)
      COMPLEX(DP), ALLOCATABLE :: dpcpart(:,:),dpcpart_J(:,:)
      TYPE(bar_type) :: barra
      INTEGER, PARAMETER :: flks(2) = [2,1]
      !
      CALL io_push_title('Compute d <a|K1d|a> / d |v>')
      !
      IF(l_nac .AND. (.NOT. PRESENT(dvg_exc_tmp_J))) &
      & CALL errore('rhs_zvector_part4','eeNAC needs state J',1)
      !
      band_group_myoffset = band_group%myoffset
      !
      ALLOCATE(z_rhs_vec_part4(npwx*npol,band_group%nlocx,kpt_pool%nloc))
      ALLOCATE(dv_vv_mat(nbndval0x-n_trunc_bands,band_group%nlocx))
      ALLOCATE(tmp_vec(npwx*npol,band_group%nlocx))
      ALLOCATE(dpcpart(npwx*npol,nbndval0x-n_trunc_bands))
      !$acc enter data create(z_rhs_vec_part4,dv_vv_mat,tmp_vec,dpcpart)
      !
      !$acc kernels present(z_rhs_vec_part4)
      z_rhs_vec_part4(:,:,:) = (0._DP,0._DP)
      !$acc end kernels
      !
      !$acc kernels present(dpcpart)
      dpcpart(:,:) = (0._DP,0._DP)
      !$acc end kernels
      !
      IF(l_nac) THEN
         !
         ALLOCATE(dpcpart_J(npwx*npol,nbndval0x-n_trunc_bands))
         !$acc enter data create(dpcpart_J)
         !
         !$acc kernels present(dpcpart_J)
         dpcpart_J(:,:) = (0._DP,0._DP)
         !$acc end kernels
         !
      ENDIF
      !
      CALL start_bar_type(barra,'zvec4',kpt_pool%nloc)
      !
      DO iks = 1,kpt_pool%nloc
         !
         IF(l_spin_flip) THEN
            iks_do = flks(iks)
         ELSE
            iks_do = iks
         ENDIF
         !
         nbndval = nbnd_occ(iks)
         flnbndval = nbnd_occ(iks_do)
         !
         nbnd_do = 0
         DO lbnd = 1,band_group%nloc
            ibnd = band_group%l2g(lbnd)+n_trunc_bands
            IF(ibnd > n_trunc_bands .AND. ibnd <= flnbndval) nbnd_do = nbnd_do+1
         ENDDO
         !
         ! ... Set k-point, spin, kinetic energy, needed by Hpsi
         !
         current_k = iks
         IF(lsda) current_spin = isk(iks)
         !
         ! ... Number of G vectors for PW expansion of wfs at k
         !
         npw = ngk(iks)
         !
         ! ... read in GS wavefunctions iks
         !
         IF(kpt_pool%nloc > 1) THEN
            IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
            CALL mp_bcast(evc,0,inter_image_comm)
            !$acc update device(evc)
         ENDIF
         !
         ! Compute the first part
         !
         IF(l_nac) THEN
            !
            ! Compute the first part
            !
            ! hybrid_kernel_term4 called once for a_I and once for a_J (derivative wrt real and complex orbitals)
            ! it uses global variable evc1_all and evc1J_all:
            ! first time evc1_all contains a_I and evc1J_all contains a_J, second time contents are switched
            !
            CALL gather_bands(dvg_exc_tmp(:,:,iks_do),evc1_all(:,:,iks_do),req)
            CALL west_mp_wait(req)
            CALL gather_bands(dvg_exc_tmp_J(:,:,iks_do),evc1J_all(:,:,iks_do),req)
            CALL west_mp_wait(req)
#if !defined(__GPU_MPI)
            !$acc update device(evc1J_all(:,:,iks_do),evc1_all(:,:,iks_do))
#endif
            !
            IF((.NOT. l_bse) .AND. l_hybrid_tddft) THEN
               CALL hybrid_kernel_term1234(current_spin,z_rhs_vec_part4(:,:,iks),l_spin_flip,4)
            ELSEIF(l_bse) THEN
               CALL bse_kernel_term4(current_spin,z_rhs_vec_part4(:,:,iks),l_spin_flip)
            ENDIF
            !
            ! switch the contents of evc1_all and evc1J_all
            !
            CALL gather_bands(dvg_exc_tmp(:,:,iks_do),evc1J_all(:,:,iks_do),req)
            CALL west_mp_wait(req)
            CALL gather_bands(dvg_exc_tmp_J(:,:,iks_do),evc1_all(:,:,iks_do),req)
            CALL west_mp_wait(req)
#if !defined(__GPU_MPI)
            !$acc update device(evc1J_all(:,:,iks_do),evc1_all(:,:,iks_do))
#endif
            !
            IF((.NOT. l_bse) .AND. l_hybrid_tddft) THEN
               CALL hybrid_kernel_term1234(current_spin,z_rhs_vec_part4(:,:,iks),l_spin_flip,4)
            ELSEIF(l_bse) THEN
               CALL bse_kernel_term4(current_spin,z_rhs_vec_part4(:,:,iks),l_spin_flip)
            ENDIF
            !
            ! the contents of evc1_all and evc1J_all are reverted back (may be unnecessary)
            !
            CALL gather_bands(dvg_exc_tmp(:,:,iks_do),evc1_all(:,:,iks_do),req)
            CALL west_mp_wait(req)
            CALL gather_bands(dvg_exc_tmp_J(:,:,iks_do),evc1J_all(:,:,iks_do),req)
            CALL west_mp_wait(req)
#if !defined(__GPU_MPI)
            !$acc update device(evc1J_all(:,:,iks_do),evc1_all(:,:,iks_do))
#endif
            !
         ELSE
            !
            IF((.NOT. l_bse) .AND. l_hybrid_tddft) THEN
               CALL hybrid_kernel_term1234(current_spin,z_rhs_vec_part4(:,:,iks),l_spin_flip,4)
            ELSEIF(l_bse) THEN
               CALL bse_kernel_term4(current_spin,z_rhs_vec_part4(:,:,iks),l_spin_flip)
            ENDIF
            !
         ENDIF
         !
         ! Compute the second part: dv_vv_mat
         !
         !$acc kernels present(tmp_vec)
         tmp_vec(:,:) = (0._DP,0._DP)
         !$acc end kernels
         !
         !$acc kernels present(dv_vv_mat)
         dv_vv_mat(:,:) = 0._DP
         !$acc end kernels
         !
         CALL bse_kernel_gamma(current_spin,evc1_all(:,:,iks),tmp_vec,l_spin_flip)
         !
         CALL glbrak_gamma(evc(:,n_trunc_bands+1:nbndval),tmp_vec,dv_vv_mat,npw,npwx,&
         & nbndval-n_trunc_bands,nbnd_do,nbndval0x-n_trunc_bands,npol)
         !
         !$acc host_data use_device(dv_vv_mat)
         CALL mp_sum(dv_vv_mat,intra_bgrp_comm)
         !$acc end host_data
         !
         IF(l_nac) THEN
            !
            !$acc host_data use_device(dvg_exc_tmp_J,dv_vv_mat,dpcpart)
            CALL DGEMM('N','T',2*npwx*npol,nbndval-n_trunc_bands,nbnd_do,-1._DP,&
            & dvg_exc_tmp_J(1,1,iks),2*npwx*npol,dv_vv_mat,nbndval0x-n_trunc_bands,0._DP,dpcpart,&
            & 2*npwx*npol)
            !$acc end host_data
            !
         ELSE
            !
            !$acc host_data use_device(dvg_exc_tmp,dv_vv_mat,dpcpart)
            CALL DGEMM('N','T',2*npwx*npol,nbndval-n_trunc_bands,nbnd_do,-1._DP,&
            & dvg_exc_tmp(1,1,iks),2*npwx*npol,dv_vv_mat,nbndval0x-n_trunc_bands,0._DP,dpcpart,&
            & 2*npwx*npol)
            !$acc end host_data
            !
         ENDIF
         !
         !$acc host_data use_device(dpcpart)
         CALL mp_sum(dpcpart,inter_bgrp_comm)
         !$acc end host_data
         !
         IF(l_nac) THEN
            !
            ! Two calls because of derivative wrt real and complex orbitals
            !
            !$acc kernels present(tmp_vec)
            tmp_vec(:,:) = (0._DP,0._DP)
            !$acc end kernels
            !
            !$acc kernels present(dv_vv_mat)
            dv_vv_mat(:,:) = 0._DP
            !$acc end kernels
            !
            CALL bse_kernel_gamma(current_spin,evc1J_all(:,:,iks),tmp_vec,l_spin_flip)
            !
            CALL glbrak_gamma(evc(:,n_trunc_bands+1:nbndval),tmp_vec,dv_vv_mat,npw,npwx,&
            & nbndval-n_trunc_bands,nbnd_do,nbndval0x-n_trunc_bands,npol)
            !
            !$acc host_data use_device(dv_vv_mat)
            CALL mp_sum(dv_vv_mat,intra_bgrp_comm)
            !$acc end host_data
            !
            !$acc host_data use_device(dvg_exc_tmp,dv_vv_mat,dpcpart_J)
            CALL DGEMM('N','T',2*npwx*npol,nbndval-n_trunc_bands,nbnd_do,-1._DP,&
            & dvg_exc_tmp(1,1,iks),2*npwx*npol,dv_vv_mat,nbndval0x-n_trunc_bands,0._DP,dpcpart_J,&
            & 2*npwx*npol)
            !$acc end host_data
            !
            !$acc host_data use_device(dpcpart_J)
            CALL mp_sum(dpcpart_J,inter_bgrp_comm)
            !$acc end host_data
            !
         ENDIF
         !
         ! compute nbnd_do for the current spin channel
         !
         nbnd_do = 0
         DO lbnd = 1,band_group%nloc
            ibnd = band_group%l2g(lbnd)+n_trunc_bands
            IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
         ENDDO
         !
         !$acc parallel loop collapse(2) present(z_rhs_vec_part4,dpcpart)
         DO lbnd = 1,nbnd_do
            DO ig = 1,npw
               !
               ! ibnd = band_group%l2g(lbnd)
               !
               ibnd = band_group_myoffset+lbnd
               !
               z_rhs_vec_part4(ig,lbnd,iks) = z_rhs_vec_part4(ig,lbnd,iks)+dpcpart(ig,ibnd)
               !
            ENDDO
         ENDDO
         !$acc end parallel
         !
         IF(l_nac) THEN
            !$acc parallel loop collapse(2) present(z_rhs_vec_part4,dpcpart_J)
            DO lbnd = 1,nbnd_do
               DO ig = 1,npw
                  !
                  ! ibnd = band_group%l2g(lbnd)
                  !
                  ibnd = band_group_myoffset+lbnd
                  !
                  z_rhs_vec_part4(ig,lbnd,iks) = z_rhs_vec_part4(ig,lbnd,iks)+dpcpart_J(ig,ibnd)
                  !
               ENDDO
            ENDDO
            !$acc end parallel
         ENDIF
         !
         IF(gstart == 2) THEN
            !$acc parallel loop present(z_rhs_vec_part4)
            DO lbnd = 1,nbnd_do
               z_rhs_vec_part4(1,lbnd,iks) = CMPLX(REAL(z_rhs_vec_part4(1,lbnd,iks),KIND=DP),KIND=DP)
            ENDDO
            !$acc end parallel
         ENDIF
         !
         !$acc parallel loop collapse(2) present(z_rhs_vec,z_rhs_vec_part4)
         DO lbnd = 1,nbnd_do
            DO ig = 1,npw
               z_rhs_vec(ig,lbnd,iks) = z_rhs_vec(ig,lbnd,iks)-z_rhs_vec_part4(ig,lbnd,iks)
            ENDDO
         ENDDO
         !$acc end parallel
         !
         CALL update_bar_type(barra,'zvec4',1)
         !
      ENDDO
      !
      CALL stop_bar_type(barra,'zvec4')
      !
      ALLOCATE(dotp(nspin))
      !
      CALL wbse_dot(z_rhs_vec_part4,z_rhs_vec_part4,band_group%nlocx,dotp)
      !
      WRITE(stdout,*)
      WRITE(stdout,"(5x,'Norm of z_rhs_vec p4 = ',ES15.8)") SUM(REAL(dotp,KIND=DP))
      !
      DEALLOCATE(dotp)
      !$acc exit data delete(z_rhs_vec_part4,dv_vv_mat,tmp_vec,dpcpart)
      DEALLOCATE(z_rhs_vec_part4)
      DEALLOCATE(dv_vv_mat)
      DEALLOCATE(tmp_vec)
      DEALLOCATE(dpcpart)
      IF(ALLOCATED(dpcpart_J)) THEN
         !$acc exit data delete(dpcpart_J)
         DEALLOCATE(dpcpart_J)
      ENDIF
      !
    END SUBROUTINE
    !
END MODULE
