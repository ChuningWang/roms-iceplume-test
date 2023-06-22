      SUBROUTINE ana_m2obc (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets 2D momentum open boundary conditions using        !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_m2obc_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     knew(ng),                                    &
     &                     GRID(ng) % angler,                           &
     &                     GRID(ng) % h,                                &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     GRID(ng) % on_u,                             &
#ifdef MASKING
     &                     GRID(ng) % umask,                            &
#endif
     &                     OCEAN(ng) % zeta)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(12)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_m2obc
!
!***********************************************************************
      SUBROUTINE ana_m2obc_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           knew,                                  &
     &                           angler, h, pm, pn, on_u,               &
#ifdef MASKING
     &                           umask,                                 &
#endif
     &                           zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: knew
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: phase

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  2D momentum open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined FJORD
!      omega=2.0_r8*pi/(12.42_r8*3600.0_r8)      ! M2 Tide period
!      omega2=2.0_r8*pi/(12.0_r8*3600.0_r8)      ! S2 Tide period
      phase=MAX(TANH(pi*(tdays(ng)-dstart)*0.1_r8), 0.0_r8)
      IF (LBC(ieast,isUbar,ng)%acquire.and.                             &
     &    LBC(ieast,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrT,JendT
          BOUNDARY(ng)%ubar_east(j)=0.0_r8
        END DO
        DO j=JstrP,JendT
          BOUNDARY(ng)%vbar_east(j)=phase*USER(1)
        END DO
      END IF

      IF (LBC(iwest,isUbar,ng)%acquire.and.                             &
     &    LBC(iwest,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrT,JendT
          BOUNDARY(ng)%ubar_west(j)=0.0_r8
        END DO
        DO j=JstrP,JendT
          BOUNDARY(ng)%vbar_west(j)=phase*USER(1)
        END DO
      END IF

      IF (LBC(isouth,isUbar,ng)%acquire.and.                            &
     &    LBC(isouth,isVbar,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrP,IendT
          BOUNDARY(ng)%ubar_south(i)=0.0_r8
        END DO
        DO i=IstrT,IendT
          BOUNDARY(ng)%vbar_south(i)=phase*USER(1)
        END DO
      END IF

      IF (LBC(inorth,isUbar,ng)%acquire.and.                            &
     &    LBC(inorth,isVbar,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrP,IendT
          BOUNDARY(ng)%ubar_north(i)=0.0_r8
        END DO
        DO i=IstrT,IendT
          BOUNDARY(ng)%vbar_north(i)=phase*USER(1)
        END DO
      END IF
#else
      IF (LBC(ieast,isUbar,ng)%acquire.and.                             &
     &    LBC(ieast,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrT,JendT
          BOUNDARY(ng)%ubar_east(j)=0.0_r8
        END DO
        DO j=JstrP,JendT
          BOUNDARY(ng)%vbar_east(j)=0.0_r8
        END DO
      END IF

      IF (LBC(iwest,isUbar,ng)%acquire.and.                             &
     &    LBC(iwest,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrT,JendT
          BOUNDARY(ng)%ubar_west(j)=0.0_r8
        END DO
        DO j=JstrP,JendT
          BOUNDARY(ng)%vbar_west(j)=0.0_r8
        END DO
      END IF

      IF (LBC(isouth,isUbar,ng)%acquire.and.                            &
     &    LBC(isouth,isVbar,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrP,IendT
          BOUNDARY(ng)%ubar_south(i)=0.0_r8
        END DO
        DO i=IstrT,IendT
          BOUNDARY(ng)%vbar_south(i)=0.0_r8
        END DO
      END IF

      IF (LBC(inorth,isUbar,ng)%acquire.and.                            &
     &    LBC(inorth,isVbar,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrP,IendT
          BOUNDARY(ng)%ubar_north(i)=0.0_r8
        END DO
        DO i=IstrT,IendT
          BOUNDARY(ng)%vbar_north(i)=0.0_r8
        END DO
      END IF
#endif

      RETURN
      END SUBROUTINE ana_m2obc_tile
