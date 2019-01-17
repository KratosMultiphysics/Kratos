/**
 *
 * @file common.h
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author David Goudin
 * @author François Pellegrini
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @author Pascal Henon
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#ifndef _spm_common_h_
#define _spm_common_h_

#include "spm.h"
#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

/********************************************************************
 * Errors functions
 */
#if defined(__GNUC__)
static inline void spm_print_error  ( const char *fmt, ...) __attribute__((format(printf,1,2)));
static inline void spm_print_warning( const char *fmt, ...) __attribute__((format(printf,1,2)));
#endif

static inline void
spm_print_error( const char *fmt, ... )
{
    va_list arglist;
    va_start(arglist, fmt);
    vfprintf(stderr, fmt, arglist);
    va_end(arglist);
}

static inline void
spm_print_warning( const char *fmt, ... )
{
    va_list arglist;
    va_start(arglist, fmt);
    fprintf(stderr, "WARNING: ");
    vfprintf(stderr, fmt, arglist);
    va_end(arglist);
}

/********************************************************************
 * CBLAS value address
 */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( a_ ) (&(a_))
#endif

/*
 * Get environment variable
 */
#if defined SPM_OS_WINDOWS

static inline int
spm_setenv( const char *var, const char *value, int overwrite ) {
    return !(SetEnvironmentVariable( var, value ));
}

static inline char *
spm_getenv( const char *var ) {
    char *str;
    int len = 512;
    int rc;
    str = (char*)malloc(len * sizeof(char));
    rc = GetEnvironmentVariable(var, str, len);
    if (rc == 0) {
        free(str);
        str = NULL;
    }
    return str;
}

static inline void
spm_cleanenv( char *str ) {
    if (str != NULL) free(str);
}

#else /* Other OS systems */

static inline int
spm_setenv( const char *var, const char *value, int overwrite ) {
    return setenv( var, value, overwrite );
}

static inline char *
spm_getenv( const char *var ) {
    return getenv( var );
}

static inline void
spm_cleanenv( char *str ) {
    (void)str;
}

#endif

#endif /* _spm_common_h_ */

