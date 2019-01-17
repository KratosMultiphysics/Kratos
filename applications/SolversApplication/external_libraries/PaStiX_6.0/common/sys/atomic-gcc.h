/*
 * Copyright (c) 2009-2015 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 */
#include <stdio.h>

static inline void pastix_mfence( void )
{
    __sync_synchronize();
}

static inline int pastix_atomic_bor_32b( volatile uint32_t* location,
                                         uint32_t value )
{
    uint32_t old_value = __sync_fetch_and_or(location, value);
    return old_value | value;
}

static inline int pastix_atomic_band_32b( volatile uint32_t* location,
                                          uint32_t value )
{
    uint32_t old_value = __sync_fetch_and_and(location, value);
    return old_value & value;
}

static inline int pastix_atomic_cas_32b( volatile uint32_t* location,
                                         uint32_t old_value,
                                         uint32_t new_value )
{
    return (__sync_bool_compare_and_swap(location, old_value, new_value) ? 1 : 0);
}

#if defined(HAVE_ATOMIC_GCC_64_BUILTINS)
static inline int pastix_atomic_cas_64b( volatile uint64_t* location,
                                         uint64_t old_value,
                                         uint64_t new_value )
{
    return (__sync_bool_compare_and_swap(location, old_value, new_value) ? 1 : 0);
}
#else
static inline int pastix_atomic_cas_64b( volatile uint64_t* location,
                                         uint64_t old_value,
                                         uint64_t new_value )
{
    fprintf(stderr, "Use of 64b CAS using atomic-gcc without __GCC_HAVE_SYNC_COMPARE_AND_SWAP_8 set\n \n");
    (void)location; (void)old_value; (void)new_value;
    return -1;
}
#endif

#if defined(HAVE_ATOMIC_GCC_128_BUILTINS)
#define PASTIX_ATOMIC_HAS_ATOMIC_CAS_128B
static inline int pastix_atomic_cas_128b( volatile __uint128_t* location,
                                         __uint128_t old_value,
                                         __uint128_t new_value )
{
    return (__sync_bool_compare_and_swap(location, old_value, new_value) ? 1 : 0);
}
#else
static inline int pastix_atomic_cas_128b( volatile uint64_t* location,
                                          uint64_t old_value,
                                          uint64_t new_value )
{
    fprintf(stderr, "Use of 128b CAS using atomic-gcc without __GCC_HAVE_SYNC_COMPARE_AND_SWAP_16 set\n \n");
    (void)location; (void)old_value; (void)new_value;
    return -1;
}
#endif

#define PASTIX_ATOMIC_HAS_ATOMIC_ADD_32B
static inline int32_t pastix_atomic_add_32b(volatile int32_t* location, int32_t i)
{
    return __sync_add_and_fetch(location, i);
}

#define PASTIX_ATOMIC_HAS_ATOMIC_SUB_32B
static inline int32_t pastix_atomic_sub_32b(volatile int32_t* location, int32_t i)
{
    return __sync_sub_and_fetch(location, i);
}
