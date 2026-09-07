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
!------------------------------------------------------------------------
SUBROUTINE do_wann()
  !------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at,alat,bg
  USE westcom,              ONLY : iuwfc,lrwfc,westpp_range,westpp_wannier_tr_rel,wannier_tr_rel,&
                                 & logfile,wann_b,wann_ng,wann_m
  USE mp_world,             ONLY : mpime,root
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : inter_image_comm,my_image_id
  USE pwcom,                ONLY : npw,current_k,ngk,nbnd
  USE control_flags,        ONLY : gamma_only
  USE buffers,              ONLY : get_buffer
  USE types_bz_grid,        ONLY : k_grid
  USE wann_loc_wfc,         ONLY : wann_init
  USE distribution_center,  ONLY : band_group
  USE class_idistribute,    ONLY : idistribute
  USE io_push,              ONLY : io_push_title
  USE json_module,          ONLY : json_file,json_core,json_value
  USE wavefunctions,        ONLY : evc
#if defined(__CUDA)
  USE west_gpu,             ONLY : allocate_gpu,deallocate_gpu
#endif
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  INTEGER :: nstate,ib
  INTEGER :: iks,il,ik
  INTEGER :: iunit
  REAL(DP) :: tmp(3)
  REAL(DP) :: wan_center(3),wan_center_cry(3)
  REAL(DP), ALLOCATABLE :: amat(:,:,:)
  REAL(DP), ALLOCATABLE :: umat(:,:)
  CHARACTER(LEN=6) :: label_k
  CHARACTER(LEN=6) :: label_b
  TYPE(json_file) :: json
  TYPE(json_core) :: jcor
  TYPE(json_value), POINTER :: jval
  !
  IF(westpp_range(2) > nbnd) CALL errore('do_wann','westpp_range(2) > nbnd',1)
  IF(.NOT. gamma_only) CALL errore('do_wann','Wannierization requires gamma_only',1)
  !
  IF(mpime == root) THEN
     CALL json%initialize()
     CALL json%load(filename=TRIM(logfile))
  ENDIF
  !
  wannier_tr_rel = westpp_wannier_tr_rel
  nstate = westpp_range(2)-westpp_range(1)+1
  !
#if defined(__CUDA)
  CALL allocate_gpu()
#endif
  !
  CALL wann_init()
  !
  ALLOCATE(amat(nstate,nstate,2*wann_ng))
  ALLOCATE(umat(nstate,nstate))
  !
  band_group = idistribute()
  CALL band_group%init(nstate,'i','westpp_range',.TRUE.)
  !
  CALL io_push_title('(B)oys/Wannier localization')
  !
  DO iks = 1,k_grid%nps ! KPOINT-SPIN LOOP
     !
     ! ... set k-point etc
     !
     current_k = iks
     npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(k_grid%nps > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
        !$acc update device(evc)
     ENDIF
     !
     CALL wbse_wann_local(westpp_range(1),westpp_range(2),amat,umat)
     !
     IF(mpime == root) THEN
        !
        WRITE(label_k,'(I6.6)') iks
        !
        ! output Wannier centers
        !
        CALL jcor%create_array(jval,'center')
        CALL json%add('output.B.K'//label_k//'.wan_center',jval)
        !
        DO ib = 1,nstate
           !
           tmp(1) = AIMAG(LOG(CMPLX(amat(ib,ib,1),amat(ib,ib,2),KIND=DP)))
           tmp(2) = AIMAG(LOG(CMPLX(amat(ib,ib,3),amat(ib,ib,4),KIND=DP)))
           tmp(3) = AIMAG(LOG(CMPLX(amat(ib,ib,5),amat(ib,ib,6),KIND=DP)))
           !
           wan_center(:) = 0._DP
           DO il = 1,3
              DO ik = 1,3
                 wan_center(il) = wan_center(il) &
                 & + tmp(ik)*wann_m(ik,il)/SQRT(wann_b(1,ik)**2+wann_b(2,ik)**2+wann_b(3,ik)**2)
              ENDDO
           ENDDO
           !
           wan_center_cry(:) = wan_center(1)*bg(1,:)/alat + wan_center(2)*bg(2,:)/alat &
           & + wan_center(3)*bg(3,:)/alat
           !
           wan_center_cry(1) = MODULO(wan_center_cry(1),1._DP)
           wan_center_cry(2) = MODULO(wan_center_cry(2),1._DP)
           wan_center_cry(3) = MODULO(wan_center_cry(3),1._DP)
           !
           wan_center(:) = wan_center_cry(1)*at(:,1)*alat + wan_center_cry(2)*at(:,2)*alat &
           & + wan_center_cry(3)*at(:,3)*alat
           !
           WRITE(label_b,'(I6)') ib
           !
           CALL json%add('output.B.K'//label_k//'.wan_center('//label_b//')',wan_center)
           !
        ENDDO
        !
        ! output transformation matrix
        !
        CALL json%add('output.B.K'//label_k//'.trans_matrix',RESHAPE(umat,[nstate*nstate]))
        !
     ENDIF
     !
  ENDDO
  !
  DEALLOCATE(amat)
  DEALLOCATE(umat)
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
  !
  IF(mpime == root) THEN
     OPEN(NEWUNIT=iunit,FILE=TRIM(logfile))
     CALL json%print(iunit)
     CLOSE(iunit)
     CALL json%destroy()
  ENDIF
  !
END SUBROUTINE
