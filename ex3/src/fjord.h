/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2017 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Estuary with Sediment Transport Test.
**
** Application flag:   ESTUARY_TEST
** Input script:       ocean_estuary_test.in
**                     sediment_estuary_test.in
*/

/* general */
#define MASKING
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS

/* iceplume */
#define ICEPLUME
#if defined ICEPLUME
# define ICEPLUME_DET_AVERAGE
# define ICEPLUME_MELT
# define ICEPLUME_MELT_TRACER
# undef ICEPLUME_WRT_AVERAGE
# undef ICEPLUME_MIX
#endif

/* tracers */
#ifdef ICEPLUME
# define T_PASSIVE
#endif

/* outputs */
#define AVERAGES
#undef PERFECT_RESTART
#undef NO_WRITE_GRID
#undef DIAGNOSTICS_UV
#undef DIAGNOSTICS_TS

/* parellel I/O */
#undef PIO_LIB

/* advection, dissipation, pressure grad, etc. */
#define UV_ADV
#define UV_COR

#define UV_VIS2
#define MIX_S_UV
#define TS_DIF2
#define MIX_GEO_TS

#ifdef SOLVE3D
# define DJ_GRADPS
#endif

/* vertical mixing */
#ifdef SOLVE3D
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#endif

#ifdef SOLVE3D
# define GLS_MIXING
# undef MY25_MIXING

# if defined GLS_MIXING || defined MY25_MIXING
#  define N2S2_HORAVG
#  define CRAIG_BANNER
#  define KANTHA_CLAYSON
#  define CHARNOK
# endif
#endif

/* horizontal mixing */
#define VISC_GRID
#define DIFF_GRID

#define ANA_DRAG
#define UV_DRAG_GRID
#define UV_LDRAG
#undef UV_QDRAG

/* tides */
#undef LTIDES
#ifdef LTIDES
# define SSH_TIDES
# define UV_TIDES
# define ADD_FSOBC
# define ADD_M2OBC
# define RAMP_TIDES
#endif

/* sediment */
#undef SEDIMENT
#ifdef SEDIMENT
# define SUSPLOAD
#endif

/* analytical functionals */
#define ANA_GRID
#define ANA_MASK
#ifndef ICEPLUME
# define ANA_INITIAL
#endif
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_M3OBC
#define ANA_TOBC

#ifdef T_PASSIVE
# define ANA_PASSIVE
#endif

#ifdef ICEPLUME
# define ANA_PSOURCE
#endif

#ifdef SEDIMENT
# define ANA_SEDIMENT
#endif

#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_SPFLUX
#define ANA_BPFLUX
