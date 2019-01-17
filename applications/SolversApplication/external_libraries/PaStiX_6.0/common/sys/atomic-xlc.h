/*
 * Copyright (c) 2012-2015 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 */
#include <stdio.h>

static inline void pastix_mfence( void )
{
    __sync();
}

static inline int pastix_atomic_bor_32b( volatile uint32_t* location,
                                         uint32_t value )
{
    uint32_t old_value = __fetch_and_or(location, value);
    return old_value | value;
}

static inline int pastix_atomic_band_32b( volatile uint32_t* location,
                                          uint32_t value )
{
    uint32_t old_value = __fetch_and_and(location, value);
    return old_value & value;
}

static inline int pastix_atomic_cas_32b( volatile uint32_t* location,
                                         uint32_t old_value,
                                         uint32_t new_value )
{
    int32_t old = (int32_t)old_value;
    return __compare_and_swap( (volatile int*)location, &old, new_value );
}

#if defined(HAVE_ATOMIC_XLC_64_BUILTINS)
static inline int pastix_atomic_cas_64b( volatile uint64_t* location,
                                         uint64_t old_value,
                                         uint64_t new_value )
{
    int64_t old = (int64_t)old_value;
    return __compare_and_swaplp( (volatile long*)location, &old, new_value );
}
#else
static inline int pastix_atomic_cas_64b( volatile uint64_t* location,
                                         uint64_t old_value,
                                         uint64_t new_value )
{
    fprintf(stderr, "Use of 64b CAS using atomic-xlc without compiler support\n \n");
    return -1;
}
#endif

#define PASTIX_ATOMIC_HAS_ATOMIC_ADD_32B
static inline int32_t pastix_atomic_add_32b( volatile int32_t *location, int32_t i )
{
    register int32_t old_val, tmp_val;

    __sync();
    do {
        old_val = __lwarx( (volatile int*)location );
        tmp_val = old_val + i;
    } while( !__stwcx( (volatile int*)location, tmp_val ) );

    return( tmp_val );
}

#define PASTIX_ATOMIC_HAS_ATOMIC_SUB_32B
static inline int32_t pastix_atomic_sub_32b( volatile int32_t *location, int32_t i )
{
    register int32_t old_val, tmp_val;

    __sync();
    do {
        old_val = __lwarx( (volatile int*)location );
        tmp_val = old_val - i;
    } while( !__stwcx( (volatile int*)location, tmp_val ) );

    return( tmp_val );
}

