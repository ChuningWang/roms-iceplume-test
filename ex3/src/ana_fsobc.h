      SUBROUTINE ana_fsobc (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets free-surface open boundary conditions using       !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_fsobc_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME( 6)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_fsobc
!
!***********************************************************************
      SUBROUTINE ana_fsobc_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS)
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
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: cff, fac, phase, val

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Free-surface open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined FJORD
!      omega=2.0_r8*pi/(12.42_r8*3600.0_r8)      ! M2 Tide period
!      omega2=2.0_r8*pi/(12.0_r8*3600.0_r8)      ! S2 Tide period
      cff=GRID(ng)%xr(Istr,Jstr)
      phase=MAX(TANH(pi*(tdays(ng)-dstart)*0.1_r8), 0.0_r8)
      IF (LBC(ieast,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrT,JendT
          fac=GRID(ng)%xr(Iend,j)-cff
          val=phase*(GRID(ng)%f(Iend,j)*fac/g)*USER(1)
          BOUNDARY(ng)%zeta_east(j)=val
        END DO
      END IF

      IF (LBC(iwest,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrT,JendT
          fac=GRID(ng)%xr(Istr,j)-cff
          val=phase*(GRID(ng)%f(Istr,j)*fac/g)*USER(1)
          BOUNDARY(ng)%zeta_west(j)=val
        END DO
      END IF

      IF (LBC(isouth,isFsur,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrT,IendT
          fac=GRID(ng)%xr(i,Jstr)-cff
          val=phase*(GRID(ng)%f(i,Jstr)*fac/g)*USER(1)
          BOUNDARY(ng)%zeta_south(i)=val
        END DO
      END IF

      IF (LBC(inorth,isFsur,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrT,IendT
          fac=GRID(ng)%xr(i,Jend)-cff
          val=phase*(GRID(ng)%f(i,Jend)*fac/g)*USER(1)
          BOUNDARY(ng)%zeta_north(i)=val
        END DO
      END IF
#else
      IF (LBC(ieast,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrT,JendT
          BOUNDARY(ng)%zeta_east(j)=0.0_r8
        END DO
      END IF

      IF (LBC(iwest,isFsur,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrT,JendT
          BOUNDARY(ng)%zeta_west(j)=0.0_r8
        END DO
      END IF

      IF (LBC(isouth,isFsur,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrT,IendT
          BOUNDARY(ng)%zeta_south(i)=0.0_r8
        END DO
      END IF

      IF (LBC(inorth,isFsur,ng)%acquire.and.                            &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrT,IendT
          BOUNDARY(ng)%zeta_north(i)=0.0_r8
        END DO
      END IF
#endif
      RETURN
      END SUBROUTINE ana_fsobc_tile
