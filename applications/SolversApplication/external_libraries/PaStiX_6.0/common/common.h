/**
 *
 * @file common.h
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author David Goudin
 * @author François Pellegrini
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @author Pascal Henon
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2018-07-16
 *
 **/
#ifndef _common_h_
#define _common_h_

#include "pastix.h"
#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include "sys/atomic.h"
#include "memory.h"
#include "integer.h"
#include "timing.h"
#include "pastixdata.h"
#include "out.h"

/*
  Macro: EXIT

  Set IPARM_ERROR_NUMBER  to module+error, dumps parameters and exit.

  Parameters:
    module - Module where the error occurs.
    error  - Value to set IPARM_ERROR_NUMBER to.
*/
#if defined(PASTIX_DEBUG_EXIT_ON_SIGSEGV)
#define EXIT(module,error) { *(int *)0 = error; }
#else
#define EXIT(module,error) { abort(); }
#endif

/********************************************************************
 * CBLAS value address
 */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( a_ ) (&(a_))
#endif

/*
 * Get environment variable
 */
#if defined(PASTIX_OS_WINDOWS)

static inline int
pastix_setenv( const char *var, const char *value, int overwrite ) {
    return !(SetEnvironmentVariable( var, value ));
}

static inline char *
pastix_getenv( const char *var ) {
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
pastix_cleanenv( char *str ) {
    if (str != NULL) free(str);
}

#else /* Other OS systems */

static inline int
pastix_setenv( const char *var, const char *value, int overwrite ) {
    return setenv( var, value, overwrite );
}

static inline char *
pastix_getenv( const char *var ) {
    return getenv( var );
}

static inline void
pastix_cleanenv( char *str ) {
    (void)str;
}

#endif


static inline int
pastix_env_is_set_to(char * str, char * value) {
    char * val;
    if ( (val = pastix_getenv(str)) &&
         !strcmp(val, value))
        return 1;
    return 0;
}

static inline int
pastix_env_is_on(char * str) {
    return pastix_env_is_set_to(str, "1");
}

static inline
int pastix_getenv_get_value_int(char * string, int default_value) {
    long int ret;
    int base = 10;
    char *endptr;
    char *str = pastix_getenv(string);
    if (str == NULL) return default_value;

    ret = strtol(str, &endptr, base);

    /* Check for various possible errors */
    if ((errno == ERANGE && (ret == LONG_MAX || ret == LONG_MIN))
        || (errno != 0 && ret == 0)) {
        perror("strtol");
        return default_value;
    }

    if (endptr == str) {
        return default_value;
    }

    if (*endptr != '\0')        /* Not necessarily an error... */
        fprintf(stderr, "Further characters after %s value: %s\n", string, endptr);
    return (int)ret;
}

/* **************************************** */

static inline void set_iparm(pastix_int_t *iparm, pastix_iparm_t offset, pastix_int_t value)
{
    if (iparm != NULL) iparm[offset] = (pastix_int_t)value;
}

static inline void set_dparm(double *dparm, pastix_dparm_t offset, double value)
{
    if (dparm != NULL) dparm[offset] = (double)value;
}

void api_dumparm(FILE *stream, pastix_int_t *iparm, double *dparm);

#endif /* _common_h_ */

