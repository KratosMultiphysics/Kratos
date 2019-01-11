/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 *
 * These macros define how C routines will be called.  ADD_ assumes that
 * they will be called by fortran, which expects C routines to have an
 * underscore postfixed to the name (Suns, and the Intel expect this).
 * NOCHANGE indicates that fortran will be calling, and that it expects
 * the name called by fortran to be identical to that compiled by the C
 * (RS6K's do this).  UPCASE says it expects C routines called by fortran
 * to be in all upcase (CRAY_PVP wants this). 
 */

#ifndef __SUPERLU_CNAMES /* allow multiple inclusions */
#define __SUPERLU_CNAMES

#define ADD_       0
#define ADD__      1
#define NOCHANGE   2
#define UPCASE     3
#define C_CALL     4

#ifdef UpCase
#define F77_CALL_C UPCASE
#endif

#ifdef NoChange
#define F77_CALL_C NOCHANGE
#endif

#ifdef Add_
#define F77_CALL_C ADD_
#endif

#ifdef Add__
#define F77_CALL_C ADD__
#endif

#ifndef F77_CALL_C
#define F77_CALL_C ADD_
#endif

#if (F77_CALL_C == ADD_)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine
 * No redefinition necessary to have following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void dgemm_(...)
 *
 * This is the default.
 */
#endif

#if (F77_CALL_C == ADD__)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine for following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void dgemm__(...)
 */
#define sgemv_    sgemv__
#define strsv_    strsv__
#define sgemm_    sgemm__
#define strsm_    strsm__

#define dgemv_    dgemv__
#define dtrsv_    dtrsv__
#define dgemm_    dgemm__
#define dtrsm_    dtrsm__

#define cgemv_    cgemv__
#define ctrsv_    ctrsv__
#define cgemm_    cgemm__
#define ctrsm_    ctrsm__

#define zgemv_    zgemv__
#define ztrsv_    ztrsv__
#define zgemm_    zgemm__
#define ztrsm_    ztrsm__

#define sasum_    sasum__
#define isamax_   isamax__
#define scopy_    scopy__
#define sscal_    sscal__
#define sger_     sger__
#define snrm2_    snrm2__
#define ssymv_    ssymv__
#define sdot_     sdot__
#define saxpy_    saxpy__
#define ssyr2_    ssyr2__
#define srot_     srot__

#define dasum_    dasum__
#define idamax_   idamax__
#define dcopy_    dcopy__
#define dscal_    dscal__
#define dgemv_    dgemv__
#define dger_     dger__
#define dnrm2_    dnrm2__
#define dsymv_    dsymv__
#define ddot_     ddot__
#define daxpy_    daxpy__
#define dsyr2_    dsyr2__
#define drot_     drot__

#define c_bridge_pdgssv_  c_bridge_pdgssv__

#endif


#if (F77_CALL_C == UPCASE)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void DGEMM(...)
 */
#define sgemv_    SGEMV
#define strsv_    STRSV
#define sgemm_    SGEMM
#define strsm_    STRSM

#define dgemv_    SGEMV
#define dtrsv_    STRSV
#define dgemm_    SGEMM
#define dtrsm_    STRSM

#define sasum_    SASUM
#define isamax_   ISAMAX
#define scopy_    SCOPY
#define sscal_    SSCAL
#define sger_     SGER
#define snrm2_    SNRM2
#define ssymv_    SSYMV
#define sdot_     SDOT
#define saxpy_    SAXPY
#define ssyr2_    SSYR2
#define srot_     SROT

#define dasum_    SASUM
#define idamax_   ISAMAX
#define dcopy_    SCOPY
#define dscal_    SSCAL
#define dgemv_    SGEMV
#define dger_     SGER
#define dnrm2_    SNRM2
#define dsymv_    SSYMV
#define ddot_     SDOT
#define daxpy_    SAXPY
#define dsyr2_    SSYR2
#define drot_     SROT

#define c_bridge_pdgssv_  C_BRIDGE_PDGSSV

#endif

#if (F77_CALL_C == NOCHANGE)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * for following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void dgemm(...)
 */
#define sgemv_    sgemv
#define strsv_    strsv
#define sgemm_    sgemm
#define strsm_    strsm

#define dgemv_    dgemv
#define dtrsv_    dtrsv
#define dgemm_    dgemm
#define dtrsm_    dtrsm

#define cgemv_    cgemv
#define ctrsv_    ctrsv
#define cgemm_    cgemm
#define ctrsm_    ctrsm

#define zgemv_    zgemv
#define ztrsv_    ztrsv
#define zgemm_    zgemm
#define ztrsm_    ztrsm

#define sasum_    sasum
#define isamax_   isamax
#define scopy_    scopy
#define sscal_    sscal
#define sger_     sger
#define snrm2_    snrm2
#define ssymv_    ssymv
#define sdot_     sdot
#define saxpy_    saxpy
#define ssyr2_    ssyr2
#define srot_     srot

#define dasum_    dasum
#define idamax_   idamax
#define dcopy_    dcopy
#define dscal_    dscal
#define dgemv_    dgemv
#define dger_     dger
#define dnrm2_    dnrm2
#define dsymv_    dsymv
#define ddot_     ddot
#define daxpy_    daxpy
#define dsyr2_    dsyr2
#define drot_     drot

#define c_bridge_pdgssv_  c_bridge_pdgssv

#endif

#endif /* __SUPERLU_CNAMES */
