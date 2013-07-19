/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2008 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

/* Atomic primitives for SiCortex machines.
 * Originally contributed by Lawrence Stewart at SiCortex.
 */

#ifndef OPA_GCC_SICORTEX_H
#define OPA_GCC_SICORTEX_H

/* FIXME do these need alignment? */
typedef struct { volatile int v;  } OPA_int_t;
typedef struct { int * volatile v; } OPA_ptr_t;

#define OPA_INT_T_INITIALIZER(val_) { (val_) }
#define OPA_PTR_T_INITIALIZER(val_) { (val_) }

#define OPA_write_barrier()      __asm__ __volatile__  ("sync" ::: "memory" )
#define OPA_read_barrier()       __asm__ __volatile__  ("sync" ::: "memory" )
#define OPA_read_write_barrier() __asm__ __volatile__  ("sync" ::: "memory" )
#define OPA_compiler_barrier()   __asm__ __volatile__  ( ""    ::: "memory" )

/* Aligned loads and stores are atomic. */
static _opa_inline int OPA_load_int(_opa_const OPA_int_t *ptr)
{
    return ptr->v;
}

/* Aligned loads and stores are atomic. */
static _opa_inline void OPA_store_int(OPA_int_t *ptr, int val)
{
    ptr->v = val;
}

/* Aligned loads and stores are atomic. */
static _opa_inline void *OPA_load_ptr(_opa_const OPA_ptr_t *ptr)
{
    return ptr->v;
}

/* Aligned loads and stores are atomic. */
static _opa_inline void OPA_store_ptr(OPA_ptr_t *ptr, void *val)
{
    ptr->v = val;
}

/* NOTE: these acquire/release operations have not been optimized, I just threw
 * down a full memory barrier.  Spending much time working on the SiCortex platform
 * doesn't really make a lot of sense since there are so few machines in
 * existence and no more will ever be built. */
static _opa_inline int OPA_load_acquire_int(_opa_const OPA_int_t *ptr)
{
    int tmp;
    tmp = ptr->v;
    OPA_read_write_barrier();
    return tmp;
}

static _opa_inline void OPA_store_release_int(OPA_int_t *ptr, int val)
{
    OPA_read_write_barrier();
    ptr->v = val;
}

static _opa_inline void *OPA_load_acquire_ptr(_opa_const OPA_ptr_t *ptr)
{
    void *tmp;
    tmp = ptr->v;
    OPA_read_write_barrier();
    return tmp;
}

static _opa_inline void OPA_store_release_ptr(OPA_ptr_t *ptr, void *val)
{
    OPA_read_write_barrier();
    ptr->v = val;
}

#include <stdint.h>

/* ICE9 rev A1 chips have a low-frequency bug that causes LL to
   fail. The workaround is to do the LL twice to make sure the data
   is in L1

   very few systems are affected

   FIXME We should either remove the workaround entirely or make it
   configurable/autodetected somehow. [goodell@ 2008/01/06]
 */
#define ICE9A_LLSC_WAR 0

/* For a description of the _opa_inline assembly constraints, see the MIPS section of:
   http://gcc.gnu.org/onlinedocs/gcc-4.3.2/gcc/Machine-Constraints.html#Machine-Constraints

   relevant excerpt:
       I - A signed 16-bit constant (for arithmetic instructions).
       J - Integer zero.

   Other _opa_inline asm knowledge worth remembering:
       r - a general register operand is allowed
       m - a memory address operand
       & - earlyclobber; operand is modified before instruction is finished using the input operands
       = - this operand is write-only

   general format:
       asm volatile ("instructions" : [outputs] : [inputs] : [clobbered_regs]);
       "memory" should be included as a clobbered reg for most of these operations.
 */

/* Atomic increment of a 32 bit value, returning the old value */
static _opa_inline int OPA_shmemi_fetch_add_4(volatile int * v, int inc)
{
        unsigned long result;
        if (ICE9A_LLSC_WAR) {
                unsigned long temp;

                __asm__ __volatile__(
                "       .set    mips3                                   \n"
                "       .set    noreorder                               \n"
                "1:     ll      %0, %2          # fetch_add_4           \n"
                "       ll      %0, %2          # fetch_add_4           \n"
                "       addu    %1, %0, %3                              \n"
                "       sc      %1, %2                                  \n"
                "       beqz    %1, 1b                                  \n"
                "        nop                                            \n"
                "       .set    reorder                                 \n"
                "       .set    mips0                                   \n"
                : "=&r" (result), "=&r" (temp), "=m" (*v)
                : "Ir" (inc)
                : "memory");
        } else {
                unsigned long temp;

                __asm__ __volatile__(
                "       .set    mips3                                   \n"
                "       .set    noreorder                               \n"
                "1:     ll      %0, %2          # fetch_add_4           \n"
                "       addu    %1, %0, %3                              \n"
                "       sc      %1, %2                                  \n"
                "       beqz    %1, 1b                                  \n"
                "        nop                                            \n"
                "       .set    reorder                                 \n"
                "       .set    mips0                                   \n"
                : "=&r" (result), "=&r" (temp), "=m" (*v)
                : "Ir" (inc)
                : "memory");
        }

        return result;
}

/* Atomic increment of a 64 bit value, returning the old value */
static _opa_inline long int OPA_shmemi_fetch_add_8(volatile long int * v, long int inc)
{
        unsigned long result;
        if (ICE9A_LLSC_WAR) {
                unsigned long temp;

                __asm__ __volatile__(
                "       .set    mips3                                   \n"
                "       .set    noreorder                               \n"
                "1:     lld     %0, %2          # fetch_add_8           \n"
                "       lld     %0, %2          # fetch_add_8           \n"
                "       daddu   %1, %0, %3                              \n"
                "       scd     %1, %2                                  \n"
                "       beqz    %1, 1b                                  \n"
                "        nop                                            \n"
                "       .set    reorder                                 \n"
                "       .set    mips0                                   \n"
                : "=&r" (result), "=&r" (temp), "=m" (*v)
                : "Ir" (inc)
                : "memory");
        } else {
                unsigned long temp;

                __asm__ __volatile__(
                "       .set    mips3                                   \n"
                "       .set    noreorder                               \n"
                "1:     lld     %0, %2          # fetch_add_8           \n"
                "       daddu   %1, %0, %3                              \n"
                "       scd     %1, %2                                  \n"
                "       beqz    %1, 1b                                  \n"
                "        nop                                            \n"
                "       .set    reorder                                 \n"
                "       .set    mips0                                   \n"
                : "=&r" (result), "=&r" (temp), "=m" (*v)
                : "Ir" (inc)
                : "memory");
        }

        return result;
}

/* Atomic swap of a 32 bit value, returning the old contents */
static _opa_inline int OPA_shmemi_swap_4(volatile int * v, int val)
{
        unsigned long result;
        if (ICE9A_LLSC_WAR) {
                unsigned long temp;

                __asm__ __volatile__(
                "       .set    mips3                                   \n"
                "       .set    noreorder                               \n"
                "1:     ll      %0, %2          # swap_4                \n"
                "       ll      %0, %2          # swap_4                \n"
                "       move    %1, %3                                  \n"
                "       sc      %1, %2                                  \n"
                "       beqz    %1, 1b                                  \n"
                "        nop                                            \n"
                "       .set    reorder                                 \n"
                "       .set    mips0                                   \n"
                : "=&r" (result), "=&r" (temp), "=m" (*v)
                : "r" (val)
                : "memory");
        } else {
                unsigned long temp;

                __asm__ __volatile__(
                "       .set    mips3                                   \n"
                "       .set    noreorder                               \n"
                "1:     ll      %0, %2          # swap_4                \n"
                "       move    %1, %3                                  \n"
                "       sc      %1, %2                                  \n"
                "       beqz    %1, 1b                                  \n"
                "        nop                                            \n"
                "       .set    reorder                                 \n"
                "       .set    mips0                                   \n"
                : "=&r" (result), "=&r" (temp), "=m" (*v)
                : "r" (val)
                : "memory");
        }

        return result;
}

/* Atomic swap of a 64 bit value, returning the old contents */
static _opa_inline long int OPA_shmemi_swap_8(volatile long int * v, long int val)
{
        unsigned long result;
        if (ICE9A_LLSC_WAR) {
                unsigned long temp;

                __asm__ __volatile__(
                "       .set    mips3                                   \n"
                "       .set    noreorder                               \n"
                "1:     lld     %0, %2          # swap_8                \n"
                "       lld     %0, %2          # swap_8                \n"
                "       move    %1, %3                                  \n"
                "       scd     %1, %2                                  \n"
                "       beqz    %1, 1b                                  \n"
                "        nop                                            \n"
                "       .set    reorder                                 \n"
                "       .set    mips0                                   \n"
                : "=&r" (result), "=&r" (temp), "=m" (*v)
                : "r" (val)
                : "memory");
        } else {
                unsigned long temp;

                __asm__ __volatile__(
                "       .set    mips3                                   \n"
                "       .set    noreorder                               \n"
                "1:     lld     %0, %2          # swap_8                \n"
                "       move    %1, %3                                  \n"
                "       scd     %1, %2                                  \n"
                "       beqz    %1, 1b                                  \n"
                "        nop                                            \n"
                "       .set    reorder                                 \n"
                "       .set    mips0                                   \n"
                : "=&r" (result), "=&r" (temp), "=m" (*v)
                : "r" (val)
                : "memory");
        }

        return result;
}

/* Atomic compare and swap of a 32 bit value, returns the old value
 * but only does the store of the val value if the old value == expect */
static _opa_inline int OPA_shmemi_cswap_4(volatile int * v, int expect, int val)
{
        unsigned long result;
        if (ICE9A_LLSC_WAR) {
                unsigned long temp;

                __asm__ __volatile__(
                "       .set    mips3                                   \n"
                "       .set    noreorder                               \n"
                "1:     ll      %0, %2          # cswap_4               \n"
                "       ll      %0, %2          # cswap_4               \n"
                "       bne     %0, %4, 1f                              \n"
                "       move    %1, %3                                  \n"
                "       sc      %1, %2                                  \n"
                "       beqz    %1, 1b                                  \n"
                "        nop                                            \n"
                "       .set    reorder                                 \n"
                "1:                                                     \n"
                "       .set    mips0                                   \n"
                : "=&r" (result), "=&r" (temp), "=m" (*v)
                : "r" (val), "Jr" (expect)
                : "memory");
        } else {
                unsigned long temp;

                __asm__ __volatile__(
                "       .set    mips3                                   \n"
                "       .set    noreorder                               \n"
                "1:     ll      %0, %2          # cswap_4               \n"
                "       bne     %0, %4, 1f                              \n"
                "       move    %1, %3                                  \n"
                "       sc      %1, %2                                  \n"
                "       beqz    %1, 1b                                  \n"
                "        nop                                            \n"
                "       .set    reorder                                 \n"
                "1:                                                     \n"
                "       .set    mips0                                   \n"
                : "=&r" (result), "=&r" (temp), "=m" (*v)
                : "r" (val), "Jr" (expect)
                : "memory");
        }

        return result;
}

/* Atomic compare and swap of a 64 bit value, returns the old value
 * but only does the store of the val value if the old value == expect */
static _opa_inline long int OPA_shmemi_cswap_8(volatile long int * v, long int expect, long int val)
{
        unsigned long result;
        if (ICE9A_LLSC_WAR) {
                unsigned long temp;

                __asm__ __volatile__(
                "       .set    mips3                                   \n"
                "       .set    noreorder                               \n"
                "1:     lld     %0, %2          # cswap_8               \n"
                "       lld     %0, %2          # cswap_8               \n"
                "       bne     %0, %4, 1f                              \n"
                "       move    %1, %3                                  \n"
                "       scd     %1, %2                                  \n"
                "       beqz    %1, 1b                                  \n"
                "        nop                                            \n"
                "       .set    reorder                                 \n"
                "1:                                                     \n"
                "       .set    mips0                                   \n"
                : "=&r" (result), "=&r" (temp), "=m" (*v)
                : "r" (val), "Jr" (expect)
                : "memory");
        } else {
                unsigned long temp;

                __asm__ __volatile__(
                "       .set    mips3                                   \n"
                "       .set    noreorder                               \n"
                "1:     lld     %0, %2          # cswap_8               \n"
                "       bne     %0, %4, 1f                              \n"
                "       move    %1, %3                                  \n"
                "       scd     %1, %2                                  \n"
                "       beqz    %1, 1b                                  \n"
                "        nop                                            \n"
                "       .set    reorder                                 \n"
                "1:                                                     \n"
                "       .set    mips0                                   \n"
                : "=&r" (result), "=&r" (temp), "=m" (*v)
                : "r" (val), "Jr" (expect)
                : "memory");
        }

        return result;
}

static _opa_inline void OPA_add_int(OPA_int_t *ptr, int val)
{
    OPA_shmemi_fetch_add_4(&ptr->v, val);
}

static _opa_inline void *OPA_cas_ptr(OPA_ptr_t *ptr, void *oldv, void *newv)
{
#if (OPA_SIZEOF_VOID_P == 8)
    return((int *) OPA_shmemi_cswap_8((volatile long int *) &ptr->v, (uintptr_t) oldv, (uintptr_t) newv));
#elif (OPA_SIZEOF_VOID_P == 4)
    return((int *) OPA_shmemi_cswap_4((volatile int *) &ptr->v, (uintptr_t) oldv, (uintptr_t) newv));
#else
#error "OPA_SIZEOF_VOID_P has an unexpected value :" OPA_QUOTE(OPA_SIZEOF_VOID_P);
#endif
}

static _opa_inline int OPA_cas_int(OPA_int_t *ptr, int oldv, int newv)
{
    return(OPA_shmemi_cswap_4(&ptr->v, oldv, newv));
}

static _opa_inline void OPA_decr_int(OPA_int_t *ptr)
{
    OPA_shmemi_fetch_add_4(&ptr->v, -1);
}

static _opa_inline int OPA_decr_and_test_int(OPA_int_t *ptr)
{
    int old = OPA_shmemi_fetch_add_4(&ptr->v, -1);
    return (old == 1);
}

static _opa_inline int OPA_fetch_and_add_int(OPA_int_t *ptr, int val)
{
    return(OPA_shmemi_fetch_add_4(&ptr->v, val));
}

static _opa_inline int OPA_fetch_and_decr_int(OPA_int_t *ptr)
{
    return(OPA_shmemi_fetch_add_4(&ptr->v, -1));
}

static _opa_inline int OPA_fetch_and_incr_int(OPA_int_t *ptr)
{
    return(OPA_shmemi_fetch_add_4(&ptr->v, 1));
}

static _opa_inline void OPA_incr_int(OPA_int_t *ptr)
{
    OPA_shmemi_fetch_add_4(&ptr->v, 1);
}

static _opa_inline int *OPA_swap_ptr(OPA_ptr_t *ptr, int *val)
{
#if (OPA_SIZEOF_VOID_P == 8)
    return((int *) OPA_shmemi_swap_8((volatile long int *) &ptr->v, (uintptr_t) val));
#elif (OPA_SIZEOF_VOID_P == 4)
    return((int *) OPA_shmemi_swap_4((volatile int *) &ptr->v, (uintptr_t) val));
#else
#error "OPA_SIZEOF_VOID_P has an unexpected value :" OPA_QUOTE(OPA_SIZEOF_VOID_P);
#endif
}

static _opa_inline int OPA_swap_int(OPA_int_t *ptr, int val)
{
    return(OPA_shmemi_swap_4(&ptr->v, val));
}

#endif /* OPA_GCC_SICORTEX_H */
