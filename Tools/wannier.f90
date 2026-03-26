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
MODULE wann_loc_wfc
  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE wann_init()
      !------------------------------------------------------------------------
      !
      USE kinds,                 ONLY : DP
      USE io_global,             ONLY : stdout
      USE constants,             ONLY : eps4,pi,tpi
      USE cell_base,             ONLY : at,alat
      USE westcom,               ONLY : wann_b,wann_g,wann_ng,wann_w,wann_m
      !
      IMPLICIT NONE
      !
      ! Workspace
      !
      REAL(DP) :: lat_a,lat_b,lat_c,ang_ab,ang_ac,ang_bc
      REAL(DP) :: det
      REAL(DP) :: m(3,3)
      INTEGER :: wann_sym
      INTEGER,PARAMETER :: SYM_UNKNOWN = 0
      INTEGER,PARAMETER :: SYM_CUBIC = 1
      INTEGER,PARAMETER :: SYM_ORTHORHOMBIC = 2
      INTEGER,PARAMETER :: SYM_HEXAGONAL = 3
      !
      wann_sym = SYM_UNKNOWN
      !
      lat_a = SQRT(SUM(at(:,1)**2))
      lat_b = SQRT(SUM(at(:,2)**2))
      lat_c = SQRT(SUM(at(:,3)**2))
      ang_ab = ACOS(DOT_PRODUCT(at(:,1),at(:,2)) / (lat_a * lat_b)) / pi * 180._DP
      ang_ac = ACOS(DOT_PRODUCT(at(:,1),at(:,3)) / (lat_a * lat_c)) / pi * 180._DP
      ang_bc = ACOS(DOT_PRODUCT(at(:,2),at(:,3)) / (lat_b * lat_c)) / pi * 180._DP
      !
      IF(ABS(ang_ab - 90._DP) < eps4 .AND. ABS(ang_ac - 90._DP) < eps4 &
      & .AND. ABS(ang_bc - 90._DP) < eps4) THEN
         IF(ABS(lat_a - lat_b) < eps4 .AND. ABS(lat_a - lat_c) < eps4 &
         & .AND. ABS(lat_b - lat_c) < eps4) THEN
            wann_sym = SYM_CUBIC
         ELSE
            wann_sym = SYM_ORTHORHOMBIC
         ENDIF
      ELSEIF(ABS(ang_ab - 120._DP) < eps4 .AND. ABS(ang_ac - 90._DP) < eps4 &
      & .AND. ABS(ang_bc - 90._DP) < eps4 .AND. ABS(lat_a - lat_b) < eps4) THEN
         wann_sym = SYM_HEXAGONAL
      ENDIF
      !
      wann_b(:,:) = 0._DP
      wann_g(:,:) = 0._DP
      wann_w(:) = 0._DP
      wann_m(:,:) = 0._DP
      !
      SELECT CASE(wann_sym)
      CASE(SYM_UNKNOWN,SYM_CUBIC)
         !
         wann_b(1,1) = tpi/alat
         wann_b(2,2) = tpi/alat
         wann_b(3,3) = tpi/alat
         wann_g(:,1) = wann_b(:,1)
         wann_g(:,2) = wann_b(:,2)
         wann_g(:,3) = wann_b(:,3)
         wann_ng = 3
         wann_w(1) = 1._DP
         wann_w(2) = 1._DP
         wann_w(3) = 1._DP
         !
         IF(wann_sym == SYM_UNKNOWN) THEN
            WRITE(stdout,"(/,7X,'** WARNING : Crystal system not implemented')")
         ELSE
            WRITE(stdout,"(/,5X,'Detected crystal system : cubic')")
         ENDIF
         !
      CASE(SYM_ORTHORHOMBIC)
         !
         wann_b(1,1) = tpi/alat
         wann_b(2,2) = tpi/alat/at(2,2)
         wann_b(3,3) = tpi/alat/at(3,3)
         wann_g(:,1) = wann_b(:,1)
         wann_g(:,2) = wann_b(:,2)
         wann_g(:,3) = wann_b(:,3)
         wann_ng = 3
         wann_w(1) = 1._DP
         wann_w(2) = at(2,2)**2
         wann_w(3) = at(3,3)**2
         !
         WRITE(stdout,"(/,5X,'Detected crystal system : orthorhombic')")
         !
      CASE(SYM_HEXAGONAL)
         !
         wann_b(1,1) = tpi/alat
         wann_b(2,1) = tpi/alat/SQRT(3._DP)
         wann_b(2,2) = tpi/alat/SQRT(3._DP)*2
         wann_b(3,3) = tpi/alat/at(3,3)
         wann_g(:,1) = wann_b(:,1)
         wann_g(:,2) = wann_b(:,2)
         wann_g(:,3) = wann_b(:,3)
         wann_g(:,4) = wann_b(:,1) - wann_b(:,2)
         wann_ng = 4
         wann_w(1) = 0.5_DP
         wann_w(2) = 0.5_DP
         wann_w(4) = 0.5_DP
         wann_w(3) = at(3,3)**2
         !
         WRITE(stdout,"(/,5X,'Detected crystal system : hexagonal')")
         !
      CASE DEFAULT
         !
         CALL errore('wann_init','unexpected wann_sym',1)
         !
      END SELECT
      !
      m(:,:) = 0._DP
      m(:,1) = wann_b(:,1) / SQRT(wann_b(1,1)**2 + wann_b(2,1)**2 + wann_b(3,1)**2)
      m(:,2) = wann_b(:,2) / SQRT(wann_b(1,2)**2 + wann_b(2,2)**2 + wann_b(3,2)**2)
      m(:,3) = wann_b(:,3) / SQRT(wann_b(1,3)**2 + wann_b(2,3)**2 + wann_b(3,3)**2)
      !
      det = m(1,1)*(m(2,2)*m(3,3) - m(2,3)*m(3,2)) &
      &   - m(1,2)*(m(2,1)*m(3,3) - m(2,3)*m(3,1)) &
      &   + m(1,3)*(m(2,1)*m(3,2) - m(2,3)*m(3,1))
      !
      wann_m(1,1) = (m(2,2)*m(3,3) - m(2,3)*m(3,2)) / det
      wann_m(1,2) = -(m(1,2)*m(3,3) - m(1,3)*m(3,2)) / det
      wann_m(1,3) = (m(1,2)*m(2,3) - m(1,3)*m(2,2)) / det
      wann_m(2,1) = -(m(2,1)*m(3,3) - m(2,3)*m(3,1)) / det
      wann_m(2,2) = (m(1,1)*m(3,3) - m(1,3)*m(3,1)) / det
      wann_m(2,3) = -(m(1,1)*m(2,3) - m(1,3)*m(2,1)) / det
      wann_m(3,1) = (m(2,1)*m(3,2) - m(2,3)*m(3,1)) / det
      wann_m(3,2) = -(m(1,1)*m(3,2) - m(1,2)*m(3,1)) / det
      wann_m(3,3) = (m(1,1)*m(2,2) - m(1,2)*m(2,1)) / det
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wann_calc_proj(proj)
      !------------------------------------------------------------------------
      !
      USE kinds,                 ONLY : DP
      USE constants,             ONLY : tpi
      USE fft_base,              ONLY : dffts
      USE scatter_mod,           ONLY : scatter_grid
      USE cell_base,             ONLY : at,alat
      USE westcom,               ONLY : wann_g,wann_ng,wann_w
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      REAL(DP),INTENT(OUT) :: proj(dffts%nnr,2*wann_ng)
      !
      ! Workspace
      !
      INTEGER :: il,ir,ix,iy,iz
      REAL(DP) :: nx,ny,nz,cry_x,cry_y,cry_z,cart_x,cart_y,cart_z
      REAL(DP) :: wc1,wc2,wc3,wc4,ws1,ws2,ws3,ws4
      !
      REAL(DP),ALLOCATABLE :: prod_gat(:)
      REAL(DP),ALLOCATABLE :: prod_distr(:)
      REAL(DP),ALLOCATABLE :: tmp(:)
      !
      proj = 0._DP
      !
      ALLOCATE(prod_gat(dffts%nr1x*dffts%nr2x*dffts%nr3x))
      ALLOCATE(prod_distr(dffts%nnr))
      ALLOCATE(tmp(4))
      !
      nx = REAL(dffts%nr1,KIND=DP)
      ny = REAL(dffts%nr2,KIND=DP)
      nz = REAL(dffts%nr3,KIND=DP)
      !
      DO il = 1,2*wann_ng
         !
         prod_gat = 0._DP
         prod_distr = 0._DP
         !
         ir = 0
         DO ix = 1,dffts%nr1
            !
            cry_x = REAL(ix-1,KIND=DP)/nx
            !
            DO iy = 1,dffts%nr2
               !
               cry_y = REAL(iy-1,KIND=DP)/ny
               !
               DO iz = 1,dffts%nr3
                  !
                  cry_z = REAL(iz-1,KIND=DP)/nz
                  !
                  cart_x = cry_x*at(1,1) + cry_y*at(1,2) + cry_z*at(1,3)
                  cart_y = cry_x*at(2,1) + cry_y*at(2,2) + cry_z*at(2,3)
                  cart_z = cry_x*at(3,1) + cry_y*at(3,2) + cry_z*at(3,3)
                  !
                  cart_x = cart_x*alat
                  cart_y = cart_y*alat
                  cart_z = cart_z*alat
                  !
                  tmp(:) = 0._DP
                  tmp(1) = wann_g(1,1)*cart_x + wann_g(2,1)*cart_y + wann_g(3,1)*cart_z
                  tmp(2) = wann_g(1,2)*cart_x + wann_g(2,2)*cart_y + wann_g(3,2)*cart_z
                  tmp(3) = wann_g(1,3)*cart_x + wann_g(2,3)*cart_y + wann_g(3,3)*cart_z
                  tmp(4) = wann_g(1,4)*cart_x + wann_g(2,4)*cart_y + wann_g(3,4)*cart_z
                  !
                  wc1 = COS(tmp(1))*SQRT(wann_w(1))
                  ws1 = SIN(tmp(1))*SQRT(wann_w(1))
                  wc2 = COS(tmp(2))*SQRT(wann_w(2))
                  ws2 = SIN(tmp(2))*SQRT(wann_w(2))
                  wc3 = COS(tmp(3))*SQRT(wann_w(3))
                  ws3 = SIN(tmp(3))*SQRT(wann_w(3))
                  wc4 = COS(tmp(4))*SQRT(wann_w(4))
                  ws4 = SIN(tmp(4))*SQRT(wann_w(4))
                  !
                  ir = (iz-1)*(dffts%nr1x*dffts%nr2x) + (iy-1)*dffts%nr1x + ix
                  IF(il == 1) prod_gat(ir) = wc1
                  IF(il == 2) prod_gat(ir) = ws1
                  IF(il == 3) prod_gat(ir) = wc2
                  IF(il == 4) prod_gat(ir) = ws2
                  IF(il == 5) prod_gat(ir) = wc3
                  IF(il == 6) prod_gat(ir) = ws3
                  IF(il == 7) prod_gat(ir) = wc4
                  IF(il == 8) prod_gat(ir) = ws4
               ENDDO
               !
            ENDDO
            !
         ENDDO
         !
         CALL scatter_grid(dffts,prod_gat,prod_distr)
         !
         DO ir = 1,dffts%nnr
            proj(ir,il) = prod_distr(ir) / (nx*ny*nz)
         ENDDO
         !
      ENDDO
      !
      DEALLOCATE(prod_gat)
      DEALLOCATE(prod_distr)
      DEALLOCATE(tmp)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wann_jade(m,a,na,u)
      !------------------------------------------------------------------------
      !
      ! Joint approximate diagonalization of eigen-matrices
      !
      ! Gygi et al., Computer Physics Communications 155, 1-6 (2003)
      !
      USE kinds,                 ONLY : DP
      USE io_global,             ONLY : stdout
      USE linear_algebra_kernel, ONLY : matdiago_dsy
      USE io_push,               ONLY : io_push_title
      USE westcom,               ONLY : wannier_tr_rel
#if defined(__CUDA)
      USE cublas
#endif
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: m
      INTEGER,INTENT(IN) :: na
      REAL(DP),INTENT(INOUT) :: a(m,m,na)
      REAL(DP),INTENT(OUT) :: u(m,m)
      !
      ! Workspace
      !
      LOGICAL :: conv
      INTEGER :: ia,i,k,mwork,iter,sweep,p,q
      REAL(DP) :: sigma,sigma_old
      REAL(DP) :: c,s,h1,h2,e1,e2,g11,g22,g12,x,y,t,tau
      !
      INTEGER,ALLOCATABLE :: top(:),bot(:)
      REAL(DP),ALLOCATABLE :: ev(:)
      REAL(DP),ALLOCATABLE :: rot(:,:),aux(:,:)
      !
      INTEGER,PARAMETER :: itermax = 100
      !
      REAL(DP) :: time_spent(2)
      REAL(DP),EXTERNAL :: get_clock
      CHARACTER(20),EXTERNAL :: human_readable_time
      !
#if defined(__CUDA)
      CALL start_clock_gpu('jade')
#else
      CALL start_clock('jade')
#endif
      !
      CALL io_push_title('Wannier (JADE)')
      !
      ALLOCATE(rot(m,m))
      ALLOCATE(aux(m,m))
      !$acc enter data create(rot,aux)
      !
      ! Handle odd m
      !
      IF(MOD(m,2) == 0) THEN
         mwork = m
      ELSE
         mwork = m+1
      ENDIF
      !
      ALLOCATE(top(mwork/2))
      ALLOCATE(bot(mwork/2))
      !
      DO k = 1,mwork/2
         top(k) = k*2 - 1
         bot(k) = k*2
      ENDDO
      !
      u(:,:) = a(:,:,1)
      !
      ALLOCATE(ev(m))
      !
      CALL matdiago_dsy(m,u,ev,.FALSE.)
      !
      DEALLOCATE(ev)
      !
      !$acc enter data copyin(a,u,top,bot)
      !
      !$acc host_data use_device(u,a,aux)
      DO ia = 1,na
         CALL DGEMM('T','N',m,m,m,1._DP,u,m,a(1,1,ia),m,0._DP,aux,m)
         CALL DGEMM('N','N',m,m,m,1._DP,aux,m,u,m,0._DP,a(1,1,ia),m)
      ENDDO
      !$acc end host_data
      !
      ! Compute initial spread
      !
      sigma_old = 0._DP
      !$acc parallel loop collapse(2) reduction(+:sigma_old) present(a) copy(sigma_old)
      DO ia = 1,na
         DO i = 1,m
            sigma_old = sigma_old + a(i,i,ia)**2
         ENDDO
      ENDDO
      !$acc end parallel
      !
      conv = .FALSE.
      !
      DO iter = 1,itermax
         !
         time_spent(1) = get_clock('jade')
         !
         DO sweep = 1,m-1
            !
            !$acc kernels present(rot)
            rot(:,:) = 0._DP
            !$acc end kernels
            !
            !$acc parallel present(top,bot,rot,a)
            !$acc loop
            DO k = 1,mwork/2
               !
               p = MIN(top(k),bot(k))
               q = MAX(top(k),bot(k))
               !
               ! Handle odd m
               !
               IF(q <= m) THEN
                  !
                  ! Compute 2x2 matrix G
                  !
                  g11 = 0._DP
                  g12 = 0._DP
                  g22 = 0._DP
                  !
                  !$acc loop seq
                  DO ia = 1,na
                     h1 = a(p,p,ia)-a(q,q,ia)
                     h2 = 2._DP*a(p,q,ia)
                     g11 = g11+h1*h1
                     g12 = g12+h1*h2
                     g22 = g22+h2*h2
                  ENDDO
                  !
                  c = 1._DP
                  s = 0._DP
                  e1 = g11
                  e2 = g22
                  !
                  ! Compute eigenvalues and eigenvectors of G
                  !
                  IF(g12*g12 > 1.E-16_DP*ABS(g11*g22)) THEN
                     tau = 0.5_DP * (g22-g11) / g12
                     t = 1.0_DP / (ABS(tau) + SQRT(1._DP+tau**2))
                     IF(tau < 0._DP) t = -t
                     c = 1._DP / SQRT(1._DP+t**2)
                     s = t * c
                     e1 = e1 - t*g12
                     e2 = e2 + t*g12
                  ENDIF
                  !
                  ! Use the eigenvector with the largest eigenvalue
                  !
                  IF(e1 > e2) THEN
                     x = c
                     y = -s
                  ELSE
                     x = s
                     y = c
                  ENDIF
                  !
                  ! Choose x >= 0 to ensure small rotation angle
                  !
                  IF(x < 0._DP) THEN
                     x = -x
                     y = -y
                  ENDIF
                  !
                  ! Compute 2x2 rotation matrix R
                  !
                  c = SQRT(0.5_DP*(x+1._DP))
                  s = y / SQRT(2._DP*(x+1._DP))
                  !
                  rot(p,p) = c
                  rot(q,p) = s
                  rot(p,q) = -s
                  rot(q,q) = c
                  !
               ELSE
                  !
                  rot(p,p) = 1._DP
                  !
               ENDIF
               !
            ENDDO
            !$acc end parallel
            !
            ! Apply rotation R to rows and columns of A
            !
            !$acc host_data use_device(rot,a,aux)
            DO ia = 1,na
               CALL DGEMM('T','N',m,m,m,1._DP,rot,m,a(1,1,ia),m,0._DP,aux,m)
               CALL DGEMM('N','N',m,m,m,1._DP,aux,m,rot,m,0._DP,a(1,1,ia),m)
            ENDDO
            !$acc end host_data
            !
            ! Accumulate unitary transformation matrix U
            !
            !$acc host_data use_device(u,rot,aux)
            CALL DGEMM('N','N',m,m,m,1._DP,u,m,rot,m,0._DP,aux,m)
            !$acc end host_data
            !
            !$acc kernels present(u,aux)
            u(:,:) = aux
            !$acc end kernels
            !
            ! Go to next round of tournament
            !
            CALL wann_tournament(top,bot,mwork)
            !
            !$acc update device(top,bot)
            !
         ENDDO
         !
         ! Compute new spread
         !
         sigma = 0._DP
         !$acc parallel loop collapse(2) reduction(+:sigma) present(a) copy(sigma)
         DO ia = 1,na
            DO i = 1,m
               sigma = sigma + a(i,i,ia)**2
            ENDDO
         ENDDO
         !$acc end parallel
         !
         WRITE(stdout,"(/,5X,'                  *----------*            *-----------------*')")
         WRITE(stdout,"(  5X,'#     Iteration = | ', I8,' |','   ','Spread = | ', ES15.8,' |')") &
         & iter, sigma
         WRITE(stdout,"(  5X,'                  *----------*            *-----------------*')")
         !
         time_spent(2) = get_clock('jade')
         !
         WRITE(stdout,"(5X,'Time spent in last iteration ',A)") &
         & TRIM(human_readable_time(time_spent(2)-time_spent(1)))
         !
         ! Check convergence
         !
         IF(ABS((sigma-sigma_old)/sigma_old) < wannier_tr_rel) THEN
            conv = .TRUE.
            EXIT
         ELSE
            sigma_old = sigma
         ENDIF
         !
      ENDDO
      !
      !$acc exit data delete(rot,aux,top,bot) copyout(a,u)
      DEALLOCATE(rot)
      DEALLOCATE(aux)
      DEALLOCATE(top)
      DEALLOCATE(bot)
      !
      IF(.NOT. conv) WRITE(stdout,'(7X,"** WARNING : JADE not converged in ",I5," steps")') itermax
      !
#if defined(__CUDA)
      CALL stop_clock_gpu('jade')
#else
      CALL stop_clock('jade')
#endif
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wann_tournament(top,bot,m)
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: m
      INTEGER,INTENT(INOUT) :: top(m/2)
      INTEGER,INTENT(INOUT) :: bot(m/2)
      !
      INTEGER,ALLOCATABLE :: new_top(:)
      INTEGER,ALLOCATABLE :: new_bot(:)
      !
      ALLOCATE(new_top(m/2))
      ALLOCATE(new_bot(m/2))
      !
      new_top(1) = top(1)
      new_top(3:m/2) = top(2:m/2-1)
      new_top(2) = bot(1)
      new_bot(1:m/2-1) = bot(2:m/2)
      new_bot(m/2) = top(m/2)
      !
      top(:) = new_top
      bot(:) = new_bot
      !
      DEALLOCATE(new_top)
      DEALLOCATE(new_bot)
      !
    END SUBROUTINE
    !
END MODULE
