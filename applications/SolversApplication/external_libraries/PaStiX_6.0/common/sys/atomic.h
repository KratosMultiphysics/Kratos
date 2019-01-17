/**
 *
 * @file atomic.h
 *
 * @copyright 2009-2012 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2014-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#ifndef _pastix_atomic_h_
#define _pastix_atomic_h_

#include <stdint.h>
#include <unistd.h>

/**
 * If the compiler provides atomic primitives we prefer to use
 * them instead of our own atomic assembly.
 */
#if defined(__FUJITSU)
  #undef HAVE_ATOMIC_XLC_32_BUILTINS
#endif
#if defined(HAVE_ATOMIC_XLC_32_BUILTINS)
#  include "common/sys/atomic-xlc.h"
#elif defined(MAC_OS_X)
#  include "common/sys/atomic-macosx.h"
#elif defined(PASTIX_ARCH_PPC)
#  if defined(__bgp__)
#    include "common/sys/atomic-ppc-bgp.h"
#  else
#    include "common/sys/atomic-ppc.h"
#  endif
#elif defined(HAVE_ATOMIC_GCC_32_BUILTINS)
#  include "common/sys/atomic-gcc.h"
#elif defined(PASTIX_ARCH_X86)
#  include "common/sys/atomic-x86_32.h"
#elif defined(PASTIX_ARCH_X86_64)
#  include "common/sys/atomic-x86_64.h"
#else
#  error "No safe atomics available"
#endif

#include <assert.h>

static inline int pastix_atomic_cas_xxb( volatile void* location,
                                         uint64_t old_value,
                                         uint64_t new_value,
                                         size_t type_size )
{
    switch(type_size){
    case 4:
        return pastix_atomic_cas_32b( (volatile uint32_t*)location,
                                     (uint32_t)old_value, (uint32_t)new_value );
    case 8:
        return pastix_atomic_cas_64b( (volatile uint64_t*)location,
                                     (uint64_t)old_value, (uint64_t)new_value );
    }
    return 0;
}

static inline uint64_t pastix_atomic_bor_xxb( volatile void* location,
                                             uint64_t or_value,
                                             size_t type_size )
{
    assert( 4 == type_size );
    (void)type_size;
    return (uint64_t)pastix_atomic_bor_32b( (volatile uint32_t*)location,
                                           (uint32_t)or_value);
}

#define pastix_atomic_band(LOCATION, OR_VALUE)  \
    (__typeof__(*(LOCATION)))pastix_atomic_band_xxb(LOCATION, OR_VALUE, sizeof(*(LOCATION)) )

#define pastix_atomic_bor(LOCATION, OR_VALUE)  \
    (__typeof__(*(LOCATION)))pastix_atomic_bor_xxb(LOCATION, OR_VALUE, sizeof(*(LOCATION)) )

#define pastix_atomic_cas(LOCATION, OLD_VALUE, NEW_VALUE)               \
    pastix_atomic_cas_xxb((volatile void*)(LOCATION),                   \
                         (uint64_t)(OLD_VALUE), (uint64_t)(NEW_VALUE), \
                         sizeof(*(LOCATION)))

#define pastix_atomic_set_mask(LOCATION, MASK) pastix_atomic_bor((LOCATION), (MASK))
#define pastix_atomic_clear_mask(LOCATION, MASK)  pastix_atomic_band((LOCATION), ~(MASK))

#ifndef PASTIX_ATOMIC_HAS_ATOMIC_INC_32B
#define PASTIX_ATOMIC_HAS_ATOMIC_INC_32B /* We now have it ! */

#ifdef PASTIX_ATOMIC_HAS_ATOMIC_ADD_32B
#define pastix_atomic_inc_32b(l)  pastix_atomic_add_32b((int32_t*)l, 1)
#else
static inline int32_t pastix_atomic_inc_32b( volatile int32_t *location )
{
    uint32_t l;
    do {
        l = (uint32_t)*location;
    } while( !pastix_atomic_cas_32b( location, l, l+1 ) );
    return (int32_t)l+1;
}
#endif  /* PASTIX_ATOMIC_HAS_ATOMIC_ADD_32B */
#endif  /* PASTIX_ATOMIC_HAS_ATOMIC_INC_32B */

#ifndef PASTIX_ATOMIC_HAS_ATOMIC_DEC_32B
#define PASTIX_ATOMIC_HAS_ATOMIC_DEC_32B /* We now have it ! */

#ifdef PASTIX_ATOMIC_HAS_ATOMIC_SUB_32B
#define pastix_atomic_dec_32b(l)  pastix_atomic_sub_32b((int32_t*)l, 1)
#else
static inline int32_t pastix_atomic_dec_32b( volatile int32_t *location )
{
    uint32_t l;
    do {
        l = (uint32_t)*location;
    } while( !pastix_atomic_cas_32b( location, l, l-1 ) );
    return (int32_t)l-1;
}
#endif  /* PASTIX_ATOMIC_HAS_ATOMIC_SUB_32B */
#endif  /* PASTIX_ATOMIC_HAS_ATOMIC_DEC_32B */

#ifndef PASTIX_ATOMIC_HAS_ATOMIC_ADD_32B
#define PASTIX_ATOMIC_HAS_ATOMIC_ADD_32B
static inline int32_t pastix_atomic_add_32b( volatile int32_t *location, int32_t d )
{
    uint32_t l, n;
    do {
        l = *location;
        n = (uint32_t)((int32_t)l + d);
    } while( !pastix_atomic_cas_32b( location, l, n ) );
    return n;
}
#endif /* PASTIX_ATOMIC_HAS_ATOMIC_ADD_32B */

typedef volatile uint32_t pastix_atomic_lock_t;
/**
 * Enumeration of lock states
 */
enum {
    PASTIX_ATOMIC_UNLOCKED = 0,
    PASTIX_ATOMIC_LOCKED   = 1
};

static inline void pastix_atomic_lock( pastix_atomic_lock_t* atomic_lock )
{
    while( !pastix_atomic_cas( atomic_lock, 0, 1) )
        /* nothing */;
}

static inline void pastix_atomic_unlock( pastix_atomic_lock_t* atomic_lock )
{
    pastix_mfence();
    *atomic_lock = 0;
}

static inline long pastix_atomic_trylock( pastix_atomic_lock_t* atomic_lock )
{
    return pastix_atomic_cas( atomic_lock, 0, 1 );
}

static inline uint64_t pastix_atomic_add_64b( volatile uint64_t *location, uint64_t d )
{
    uint64_t l, n;
    do {
        l = (*location);
        n = (l + d);
    } while( !pastix_atomic_cas_64b( location, l, n ) );
    return n;
}
#endif /* _pastix_atomic_h_ */
