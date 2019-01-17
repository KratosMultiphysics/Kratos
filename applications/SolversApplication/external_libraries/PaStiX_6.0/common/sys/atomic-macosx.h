/*
 * Copyright (c) 2009-2014 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 */

#ifndef MAC_OS_X
#error This file should only be included on MAC OS X (Snow Leopard
#endif

#include <libkern/OSAtomic.h>

static inline void pastix_mfence( void )
{
    OSMemoryBarrier();
}

static inline int pastix_atomic_bor_32b( volatile uint32_t* location,
                                          uint32_t value )
{
    return OSAtomicOr32( value, location );
}

static inline int pastix_atomic_band_32b( volatile uint32_t* location,
                                          uint32_t value )
{
    return OSAtomicAnd32( value, location );
}

static inline int pastix_atomic_cas_32b( volatile uint32_t* location,
                                          uint32_t old_value,
                                          uint32_t new_value )
{
    return OSAtomicCompareAndSwap32( old_value, new_value, (volatile int32_t*)location );
}

static inline int pastix_atomic_cas_64b( volatile uint64_t* location,
                                          uint64_t old_value,
                                          uint64_t new_value )
{
    return OSAtomicCompareAndSwap64( old_value, new_value, (volatile int64_t*)location );
}

#define PASTIX_ATOMIC_HAS_ATOMIC_INC_32B
static inline int32_t pastix_atomic_inc_32b( volatile int32_t *location )
{
    return OSAtomicIncrement32( (int32_t*)location );
}

#define PASTIX_ATOMIC_HAS_ATOMIC_DEC_32B
static inline int32_t pastix_atomic_dec_32b( volatile int32_t *location )
{
    return OSAtomicDecrement32( (int32_t*)location );
}

#define PASTIX_ATOMIC_HAS_ATOMIC_ADD_32B
static inline int32_t pastix_atomic_add_32b( volatile int32_t *location, int32_t i )
{
    return OSAtomicAdd32( i, location );
}

#define PASTIX_ATOMIC_HAS_ATOMIC_SUB_32B
static inline int32_t pastix_atomic_sub_32b( volatile int32_t *location, int32_t i )
{
    return OSAtomicAdd32( -i, location );
}

#if defined(HAVE_ATOMIC_GCC_128_BUILTINS)
#define PASTIX_ATOMIC_HAS_ATOMIC_CAS_128B
static inline int pastix_atomic_cas_128b( volatile __int128_t* location,
                                         __int128_t old_value,
                                         __int128_t new_value )
{
    return (__sync_bool_compare_and_swap(location, old_value, new_value) ? 1 : 0);
}
#endif
