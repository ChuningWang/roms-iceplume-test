      SUBROUTINE ana_grid (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets model grid using an analytical expressions.       !
!                                                                      !
!  On Output:  stored in common blocks:                                !
!                                                                      !
!                           "grid"    (file grid.h)                    !
!                           "scalars" (file scalar.h)                  !
!                                                                      !
!     el       Length (m) of domain box in the ETA-direction.          !
!     f        Coriolis parameter (1/seconds) at RHO-points.           !
!     h        Bathymetry (meters; positive) at RHO-points.            !
!     hmin     Minimum depth of bathymetry (m).                        !
!     hmax     Maximum depth of bathymetry (m).                        !
!     pm       Coordinate transformation metric "m" (1/meters)         !
!              associated with the differential distances in XI        !
!              at RHO-points.                                          !
!     pn       Coordinate transformation metric "n" (1/meters)         !
!              associated with the differential distances in ETA.      !
!              at RHO-points.                                          !
!     xl       Length (m) of domain box in the XI-direction.           !
!     xp       XI-coordinates (m) at PSI-points.                       !
!     xr       XI-coordinates (m) at RHO-points.                       !
!     yp       ETA-coordinates (m) at PSI-points.                      !
!     yr       ETA-coordinates (m) at RHO-points.                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_grid_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    GRID(ng) % angler,                            &
#if defined CURVGRID && defined UV_ADV
     &                    GRID(ng) % dmde,                              &
     &                    GRID(ng) % dndx,                              &
#endif
#ifdef ICESHELF
     &                    GRID(ng) % zice,                              &
#endif
#ifdef SPHERICAL
     &                    GRID(ng) % lonp,                              &
     &                    GRID(ng) % lonr,                              &
     &                    GRID(ng) % lonu,                              &
     &                    GRID(ng) % lonv,                              &
     &                    GRID(ng) % latp,                              &
     &                    GRID(ng) % latr,                              &
     &                    GRID(ng) % latu,                              &
     &                    GRID(ng) % latv,                              &
#else
     &                    GRID(ng) % xp,                                &
     &                    GRID(ng) % xr,                                &
     &                    GRID(ng) % xu,                                &
     &                    GRID(ng) % xv,                                &
     &                    GRID(ng) % yp,                                &
     &                    GRID(ng) % yr,                                &
     &                    GRID(ng) % yu,                                &
     &                    GRID(ng) % yv,                                &
#endif
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % f,                                 &
     &                    GRID(ng) % h)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME( 7)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_grid
!
!***********************************************************************
      SUBROUTINE ana_grid_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          angler,                                 &
#if defined CURVGRID && defined UV_ADV
     &                          dmde, dndx,                             &
#endif
#ifdef ICESHELF
     &                          zice,                                   &
#endif
#ifdef SPHERICAL
     &                          lonp, lonr, lonu, lonv,                 &
     &                          latp, latr, latu, latv,                 &
#else
     &                          xp, xr, xu, xv,                         &
     &                          yp, yr, yu, yv,                         &
#endif
     &                          pn, pm, f, h)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_iounits
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
      USE stats_mod, ONLY : stats_2dfld
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: angler(LBi:,LBj:)
# if defined CURVGRID && defined UV_ADV
      real(r8), intent(out) :: dmde(LBi:,LBj:)
      real(r8), intent(out) :: dndx(LBi:,LBj:)
# endif
# ifdef ICESHELF
      real(r8), intent(out) :: zice(LBi:,LBj:)
# endif
# ifdef SPHERICAL
      real(r8), intent(out) :: lonp(LBi:,LBj:)
      real(r8), intent(out) :: lonr(LBi:,LBj:)
      real(r8), intent(out) :: lonu(LBi:,LBj:)
      real(r8), intent(out) :: lonv(LBi:,LBj:)
      real(r8), intent(out) :: latp(LBi:,LBj:)
      real(r8), intent(out) :: latr(LBi:,LBj:)
      real(r8), intent(out) :: latu(LBi:,LBj:)
      real(r8), intent(out) :: latv(LBi:,LBj:)
# else
      real(r8), intent(out) :: xp(LBi:,LBj:)
      real(r8), intent(out) :: xr(LBi:,LBj:)
      real(r8), intent(out) :: xu(LBi:,LBj:)
      real(r8), intent(out) :: xv(LBi:,LBj:)
      real(r8), intent(out) :: yp(LBi:,LBj:)
      real(r8), intent(out) :: yr(LBi:,LBj:)
      real(r8), intent(out) :: yu(LBi:,LBj:)
      real(r8), intent(out) :: yv(LBi:,LBj:)
# endif
      real(r8), intent(out) :: pn(LBi:,LBj:)
      real(r8), intent(out) :: pm(LBi:,LBj:)
      real(r8), intent(out) :: f(LBi:,LBj:)
      real(r8), intent(out) :: h(LBi:,LBj:)
#else
      real(r8), intent(out) :: angler(LBi:UBi,LBj:UBj)
# if defined CURVGRID && defined UV_ADV
      real(r8), intent(out) :: dmde(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: dndx(LBi:UBi,LBj:UBj)
# endif
# ifdef ICESHELF
      real(r8), intent(out) :: zice(LBi:UBi,LBj:UBj)
# endif
# ifdef SPHERICAL
      real(r8), intent(out) :: lonp(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: lonu(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: lonv(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: latp(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: latr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: latu(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: latv(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(out) :: xp(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: xu(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: xv(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: yp(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: yr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: yu(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: yv(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: f(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: h(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      logical, save :: first = .TRUE.

      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, ival, j, k, icen, iwid

      real(r8), parameter :: twopi = 2.0_r8*pi

      real(r8) :: Esize, Xsize, beta, cff, depth, dth
      real(r8) :: dx, dy, f0, r, theta, val1, val2

      real(r8) :: wrkX(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: wrkY(IminS:ImaxS,JminS:JmaxS)

      TYPE (T_STATS), save :: Stats(16)

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set grid parameters:
!
!     Xsize    Length (m) of domain box in the XI-direction.
!     Esize    Length (m) of domain box in the ETA-direction.
!     depth    Maximum depth of bathymetry (m).
!     f0       Coriolis parameter, f-plane constant (1/s).
!     beta     Coriolis parameter, beta-plane constant (1/s/m).
!-----------------------------------------------------------------------
!
#if defined FJORD
      icen=NINT(REAL(Mm(ng)+1,r8)/2.0_r8)
      iwid=NINT(REAL(Mm(ng)-48-1,r8)/2.0_r8)
      Xsize=0.0_r8
      DO i=1,Lm(ng)
        IF (i.le.201) THEN
          dx=300.0_r8
        ELSE
          dx=300.0_r8+50.0_r8*REAL((i-201),r8)
        ENDIF
        Xsize=Xsize+dx
      ENDDO
      Esize=0.0_r8
      DO j=1,Mm(ng)
        ival=ABS(j-icen)
        IF (ival.le.iwid) THEN
          dy=220.0_r8
        ELSE
          dy=220.0_r8+50.0_r8*REAL(ival-iwid,r8)
        ENDIF
        Esize=Esize+dy
      ENDDO
      depth=400.0_r8
      f0=USER(2)*1.379e-4_r8
      beta=0.0_r8
#else
      ana_grid.h: no values provided for Xsize, Esize, depth, f0, beta.
#endif
!
!  Load grid parameters to global storage.
!
      IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
        xl(ng)=Xsize
        el(ng)=Esize
      END IF
!
!-----------------------------------------------------------------------
!  Initialize field statistics structure.
!-----------------------------------------------------------------------
!
      IF (first) THEN
        first=.FALSE.
        DO i=1,SIZE(Stats,1)
          Stats(i) % count=0.0_r8
          Stats(i) % min=Large
          Stats(i) % max=-Large
          Stats(i) % avg=0.0_r8
          Stats(i) % rms=0.0_r8
        END DO
      END IF
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) WRITE (stdout,'(1x)')
!
!-----------------------------------------------------------------------
!  Compute the (XI,ETA) coordinates at PSI- and RHO-points.
!  Set grid spacing (m).
!-----------------------------------------------------------------------
!
!  Determine I- and J-ranges for computing grid data.  These ranges
!  are special in periodic boundary conditons since periodicity cannot
!  be imposed in the grid coordinates.
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=Istr-1
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=Iend+1
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=Jstr-1
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=Jend+1
      ELSE
        Jmax=Jend
      END IF
!
      icen=NINT(REAL(Mm(ng)+1,r8)/2.0_r8)
      iwid=NINT(REAL(Mm(ng)-48-1,r8)/2.0_r8)
      DO j=Jmin,Jmax
        DO i=Imin,Imax
!
          ival=j-icen
          IF (i.le.201) THEN
            dx=300.0_r8
          ELSE
            dx=300.0_r8+50.0_r8*REAL((i-201),r8)
          ENDIF
          IF (ABS(ival).le.iwid) THEN
            dy=220.0_r8
          ELSE
            dy=220.0_r8+50.0_r8*REAL(ABS(ival)-iwid,r8)
          ENDIF
!
          IF (i.le.201) THEN
            xr(i,j)=300.0_r8*REAL(i-1,r8)-150.0_r8
          ELSE
            xr(i,j)=300.0_r8*REAL(i-1,r8)+50.0_r8*                      &
     &        (REAL(i-201,r8)+REAL((i-201)*(i-202),r8)/2.0_r8)-         &
     &        25.0_r8*REAL(i-201)-150.0_r8
          ENDIF
          IF (ABS(ival).le.iwid) THEN
            yr(i,j)=220.0_r8*REAL(ival)
          ELSE IF (ival.gt.iwid) THEN
            yr(i,j)=220.0_r8*REAL(ival)+50.0_r8*                        &
     & (REAL(ival-iwid,r8)+                                             &
     &  REAL((ival-iwid)*(ival-iwid-1),r8)/2.0_r8)-                     &
     & 25.0_r8*REAL(ival-iwid)
          ELSE IF (ival.lt.-1*iwid) THEN
            yr(i,j)=220.0_r8*REAL(ival)-50.0_r8*                        &
     & (-1*REAL(ival+iwid,r8)+                                          &
     &  REAL((ival+iwid)*(ival+iwid+1),r8)/2.0_r8)-                     &
     & 25.0_r8*REAL(ival+iwid)
          ENDIF
          xp(i,j)=xr(i,j)-0.5_r8*dx
          xu(i,j)=xp(i,j)
          xv(i,j)=xr(i,j)
          yp(i,j)=yr(i,j)-0.5_r8*dy
          yu(i,j)=yr(i,j)
          yv(i,j)=yp(i,j)
        END DO
      END DO
!
!  Report statistics.
!
#ifdef SPHERICAL
      CALL stats_2dfld (ng, tile, iNLM, p2dvar, Stats(1),               &
     &                  LBi, UBi, LBj, UBj, lonp)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'longitude of PSI-points: lon_psi',           &
     &                     ng, Stats(1)%min, Stats(1)%max
      END IF
      CALL stats_2dfld (ng, tile, iNLM, p2dvar, Stats(2),               &
     &                  LBi, UBi, LBj, UBj, latp)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'latitude of PSI-points: lat_psi',            &
     &                     ng, Stats(2)%min, Stats(2)%max
      END IF

      CALL stats_2dfld (ng, tile, iNLM, r2dvar, Stats(3),               &
     &                  LBi, UBi, LBj, UBj, lonr)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'longitude of RHO-points: lon_rho',           &
     &                     ng, Stats(3)%min, Stats(3)%max
      END IF
      CALL stats_2dfld (ng, tile, iNLM, r2dvar, Stats(4),               &
     &                  LBi, UBi, LBj, UBj, latr)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'latitude of RHO-points: lat_rho',            &
     &                     ng, Stats(4)%min, Stats(4)%max
      END IF

      CALL stats_2dfld (ng, tile, iNLM, u2dvar, Stats(5),               &
     &                  LBi, UBi, LBj, UBj, lonu)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'longitude of U-points: lon_u',               &
     &                     ng, Stats(5)%min, Stats(5)%max
      END IF
      CALL stats_2dfld (ng, tile, iNLM, u2dvar, Stats(6),               &
     &                  LBi, UBi, LBj, UBj, latu)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'latitude of U-points: lat_u',                &
     &                     ng, Stats(6)%min, Stats(6)%max
      END IF

      CALL stats_2dfld (ng, tile, iNLM, v2dvar, Stats(7),               &
     &                  LBi, UBi, LBj, UBj, lonv)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'longitude of V-points: lon_v',               &
     &                     ng, Stats(7)%min, Stats(7)%max
      END IF
      CALL stats_2dfld (ng, tile, iNLM, v2dvar, Stats(8),               &
     &                  LBi, UBi, LBj, UBj, latv)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'latitude of V-points: lat_v',                &
     &                     ng, Stats(8)%min, Stats(8)%max
      END IF
#else
      CALL stats_2dfld (ng, tile, iNLM, p2dvar, Stats(1),               &
     &                  LBi, UBi, LBj, UBj, xp)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'x-location of PSI-points: x_psi',            &
     &                     ng, Stats(1)%min, Stats(1)%max
      END IF
      CALL stats_2dfld (ng, tile, iNLM, p2dvar, Stats(2),               &
     &                  LBi, UBi, LBj, UBj, yp)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'y-location of PSI-points: y_psi',            &
     &                     ng, Stats(2)%min, Stats(2)%max
      END IF

      CALL stats_2dfld (ng, tile, iNLM, r2dvar, Stats(3),               &
     &                  LBi, UBi, LBj, UBj, xr)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'x-location of RHO-points: x_rho',            &
     &                     ng, Stats(3)%min, Stats(3)%max
      END IF
      CALL stats_2dfld (ng, tile, iNLM, r2dvar, Stats(4),               &
     &                  LBi, UBi, LBj, UBj, yr)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'y-location of RHO-points: y_rho',            &
     &                     ng, Stats(4)%min, Stats(4)%max
      END IF

      CALL stats_2dfld (ng, tile, iNLM, u2dvar, Stats(5),               &
     &                  LBi, UBi, LBj, UBj, xu)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'x-location of U-points: x_u',                &
     &                     ng, Stats(5)%min, Stats(5)%max
      END IF
      CALL stats_2dfld (ng, tile, iNLM, u2dvar, Stats(6),               &
     &                  LBi, UBi, LBj, UBj, yu)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'y-location of U-points: y_u',                &
     &                     ng, Stats(6)%min, Stats(6)%max
      END IF

      CALL stats_2dfld (ng, tile, iNLM, v2dvar, Stats(7),               &
     &                  LBi, UBi, LBj, UBj, xv)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'x-location of V-points: x_v',                &
     &                     ng, Stats(7)%min, Stats(7)%max
      END IF
      CALL stats_2dfld (ng, tile, iNLM, v2dvar, Stats(8),               &
     &                  LBi, UBi, LBj, UBj, yv)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'y-location of V-points: y_v',                &
     &                     ng, Stats(8)%min, Stats(8)%max
      END IF
#endif

#ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
# ifdef SPHERICAL
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    lonp, lonr, lonu, lonv)
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    latp, latr, latu, latv)
# else
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    xp, xr, xu, xv)
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    yp, yr, yu, yv)
# endif
#endif
!
!-----------------------------------------------------------------------
! Compute coordinate transformation metrics at RHO-points "pm" and
! "pn"  (1/m) associated with the differential distances in XI and
! ETA, respectively.
!-----------------------------------------------------------------------
!
#define J_RANGE MIN(JstrT,Jstr-1),MAX(Jend+1,JendT)
#define I_RANGE MIN(IstrT,Istr-1),MAX(Iend+1,IendT)

      icen=NINT(REAL(Mm(ng)+1,r8)/2.0_r8)
      iwid=NINT(REAL(Mm(ng)-48-1,r8)/2.0_r8)
      DO j=J_RANGE
        DO i=I_RANGE
!
          ival=ABS(j-icen)
          IF (i.le.201) THEN
            dx=300.0_r8
          ELSE
            dx=300.0_r8+50.0_r8*REAL((i-201),r8)
          ENDIF
          IF (ival.le.iwid) THEN
            dy=220.0_r8
          ELSE
            dy=220.0_r8+50.0_r8*REAL(ival-iwid,r8)
          ENDIF
!
          wrkX(i,j)=1.0_r8/dx
          wrkY(i,j)=1.0_r8/dy
        END DO
      END DO
#undef J_RANGE
#undef I_RANGE
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          pm(i,j)=wrkX(i,j)
          pn(i,j)=wrkY(i,j)
        END DO
      END DO
!
!  Report statistics.
!
      CALL stats_2dfld (ng, tile, iNLM, r2dvar, Stats(9),               &
     &                  LBi, UBi, LBj, UBj, pm)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'reciprocal XI-grid spacing: pm',             &
     &                     ng, Stats(9)%min, Stats(9)%max
      END IF
      CALL stats_2dfld (ng, tile, iNLM, r2dvar, Stats(10),              &
     &                  LBi, UBi, LBj, UBj, pn)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'reciprocal ETA-grid spacing: pn',            &
     &                     ng, Stats(10)%min, Stats(10)%max
      END IF
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pm)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pn)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pm, pn)
#endif

#if (defined CURVGRID && defined UV_ADV)
!
!-----------------------------------------------------------------------
!  Compute d(1/n)/d(xi) and d(1/m)/d(eta) at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          dndx(i,j)=0.5_r8*((1.0_r8/wrkY(i+1,j  ))-                     &
     &                      (1.0_r8/wrkY(i-1,j  )))
          dmde(i,j)=0.5_r8*((1.0_r8/wrkX(i  ,j+1))-                     &
     &                      (1.0_r8/wrkX(i  ,j-1)))
        END DO
      END DO
!
!  Report statistics.
!
      CALL stats_2dfld (ng, tile, iNLM, r2dvar, Stats(11),              &
     &                  LBi, UBi, LBj, UBj, dmde)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'ETA-derivative of inverse metric '//         &
     &                    'factor pm: dmde',                            &
     &                     ng, Stats(11)%min, Stats(11)%max
      END IF
      CALL stats_2dfld (ng, tile, iNLM, r2dvar, Stats(12),              &
     &                  LBi, UBi, LBj, UBj, dndx)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'XI-derivative of inverse metric '//          &
     &                    'factor pn: dndx',                            &
     &                     ng, Stats(12)%min, Stats(12)%max
      END IF
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          dndx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          dmde)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    dndx, dmde)
# endif
#endif
!
!-----------------------------------------------------------------------
! Angle (radians) between XI-axis and true EAST at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          angler(i,j)=0.0_r8
        END DO
      END DO
!
!  Report Statistics.
!
      CALL stats_2dfld (ng, tile, iNLM, r2dvar, Stats(13),              &
     &                  LBi, UBi, LBj, UBj, angler)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'angle between XI-axis and EAST: '//          &
     &                    'angler',                                     &
     &                     ng, Stats(13)%min, Stats(13)%max
      END IF
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          angler)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    angler)
#endif
!
!-----------------------------------------------------------------------
!  Compute Coriolis parameter (1/s) at RHO-points.
!-----------------------------------------------------------------------
!
      val1=0.5_r8*Esize
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          f(i,j)=f0+beta*(yr(i,j)-val1)
        END DO
      END DO
!
!  Report Statistics.
!
      CALL stats_2dfld (ng, tile, iNLM, r2dvar, Stats(14),              &
     &                  LBi, UBi, LBj, UBj, f)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'Coriolis parameter at RHO-points: f',        &
     &                     ng, Stats(14)%min, Stats(14)%max
      END IF
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          f)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    f)
#endif
!
!-----------------------------------------------------------------------
!  Set bathymetry (meters; positive) at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          h(i,j)=depth
        END DO
      END DO
!
!  Report Statistics.
!
      CALL stats_2dfld (ng, tile, iNLM, r2dvar, Stats(15),              &
     &                  LBi, UBi, LBj, UBj, h)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'bathymetry at RHO-points: h',                &
     &                     ng, Stats(15)%min, Stats(15)%max
      END IF
      hmin(ng)=Stats(15)%min
      hmax(ng)=Stats(15)%max
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          h)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    h)
#endif
#ifdef ICESHELF
!
!-----------------------------------------------------------------------
!  Set depth of ice shelf (meters; negative) at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          zice(i,j)=0.0_r8
        END DO
      END DO
!
!  Report Statistics.
!
      CALL stats_2dfld (ng, tile, iNLM, r2dvar, Stats(16),              &
     &                  LBi, UBi, LBj, UBj, zice)
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        WRITE (stdout,10) 'ice shelf thickness: zice',                  &
     &                     ng, Stats(16)%min, Stats(16)%max
      END IF
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zice)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zice)
# endif
#endif
!
  10  FORMAT (3x,' ANA_GRID    - ',a,/,19x,                             &
     &        '(Grid = ',i2.2,', Min = ',1p,e15.8,0p,                   &
     &                         ' Max = ',1p,e15.8,0p,')')

      RETURN
      END SUBROUTINE ana_grid_tile
