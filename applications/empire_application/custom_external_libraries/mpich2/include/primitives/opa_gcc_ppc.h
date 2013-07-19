/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*  
 *  (C) 2008 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef OPA_GCC_PPC_H_INCLUDED
#define OPA_GCC_PPC_H_INCLUDED

/* these need to be aligned on an 8-byte boundary to work on a BG/P */
typedef struct { volatile int v    OPA_ATTRIBUTE((aligned (8))); } OPA_int_t;
typedef struct { void * volatile v OPA_ATTRIBUTE((aligned (8))); } OPA_ptr_t;

#define OPA_INT_T_INITIALIZER(val_) { (val_) }
#define OPA_PTR_T_INITIALIZER(val_) { (val_) }

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

/* a useful link for PPC memory ordering issues:
 *   http://www.rdrop.com/users/paulmck/scalability/paper/N2745r.2009.02.22a.html
 *
 * lwsync: orders L-L, S-S, L-S, but *not* S-L (i.e. gives x86-ish ordering)
 * eieio: orders S-S (but only for cacheable memory, not for MMIO)
 * sync: totally orders memops
 * isync: force all preceeding insns to appear complete before starting
 *        subsequent insns, but w/o cumulativity (very confusing)
 */
/* helper macros */
#define OPA_ppc_lwsync_() __asm__ __volatile__  ( "lwsync" ::: "memory" )
#define OPA_ppc_hwsync_() __asm__ __volatile__  ( "sync"   ::: "memory" )
#define OPA_ppc_eieio_()  __asm__ __volatile__  ( "eieio"  ::: "memory" )

#define OPA_write_barrier()      OPA_ppc_eieio_()
#define OPA_read_barrier()       OPA_ppc_lwsync_()
#define OPA_read_write_barrier() OPA_ppc_hwsync_()
#define OPA_compiler_barrier()   __asm__ __volatile__  ( "" ::: "memory" )

/* NOTE-PPC-1 we use lwsync, although I think we might be able to use
 * conditional-branch+isync in some cases (load_acquire?) once we understand it
 * better */

static _opa_inline int   OPA_load_acquire_int(_opa_const OPA_int_t *ptr)
{
    int tmp;
    tmp = ptr->v;
    OPA_ppc_lwsync_(); /* NOTE-PPC-1 */
    return tmp;
}
static _opa_inline void  OPA_store_release_int(OPA_int_t *ptr, int val)
{
    OPA_ppc_lwsync_();
    ptr->v = val;
}
static _opa_inline void *OPA_load_acquire_ptr(_opa_const OPA_ptr_t *ptr)
{
    void *tmp;
    tmp = ptr->v;
    OPA_ppc_lwsync_(); /* NOTE-PPC-1 */
    return tmp;
}
static _opa_inline void  OPA_store_release_ptr(OPA_ptr_t *ptr, void *val)
{
    OPA_ppc_lwsync_();
    ptr->v = val;
}


/*
   load-link/store-conditional (LL/SC) primitives.  We LL/SC implement
   these here, which are arch-specific, then use the generic
   implementations from opa_emulated.h */

static _opa_inline int OPA_LL_int(OPA_int_t *ptr)
{
    int val;
    __asm__ __volatile__ ("lwarx %[val],0,%[ptr]"
                          : [val] "=r" (val)
                          : [ptr] "r" (&ptr->v)
                          : "cc");

    return val;
}

/* Returns non-zero if the store was successful, zero otherwise. */
static _opa_inline int OPA_SC_int(OPA_int_t *ptr, int val)
{
    int ret = 1; /* init to non-zero, will be reset to 0 if SC was successful */
    __asm__ __volatile__ ("stwcx. %[val],0,%[ptr];\n"
                          "beq 1f;\n"
                          "li %[ret], 0;\n"
                          "1: ;\n"
                          : [ret] "=r" (ret)
                          : [ptr] "r" (&ptr->v), [val] "r" (val), "0" (ret)
                          : "cc", "memory");
    return ret;
}


/* Pointer versions of LL/SC. */

/* Set OPA_SS (Size Suffix) which is used to choose between lwarx/stwcx and
 * ldarx/stdcx when using 4 or 8 byte pointer operands */
#if OPA_SIZEOF_VOID_P == 4
#define OPA_SS "w"
#elif OPA_SIZEOF_VOID_P == 8
#define OPA_SS "d"
#else
#error OPA_SIZEOF_VOID_P is not 4 or 8
#endif


static _opa_inline void *OPA_LL_ptr(OPA_ptr_t *ptr)
{
    void *val;
    __asm__ __volatile__ ("l"OPA_SS"arx %[val],0,%[ptr]"
                          : [val] "=r" (val)
                          : [ptr] "r" (&ptr->v)
                          : "cc");

    return val;
}

/* Returns non-zero if the store was successful, zero otherwise. */
static _opa_inline int OPA_SC_ptr(OPA_ptr_t *ptr, void *val)
{
    int ret = 1; /* init to non-zero, will be reset to 0 if SC was successful */
    __asm__ __volatile__ ("st"OPA_SS"cx. %[val],0,%[ptr];\n"
                          "beq 1f;\n"
                          "li %[ret], 0;\n"
                          "1: ;\n"
                          : [ret] "=r" (ret)
                          : [ptr] "r" (&ptr->v), [val] "r" (val), "0" (ret)
                          : "cc", "memory");
    return ret;
}

#undef OPA_SS

/* necessary to enable LL/SC emulation support */
#define OPA_LL_SC_SUPPORTED 1

/* Implement all function using LL/SC */
#define OPA_add_int_by_llsc            OPA_add_int
#define OPA_incr_int_by_llsc           OPA_incr_int
#define OPA_decr_int_by_llsc           OPA_decr_int
#define OPA_decr_and_test_int_by_llsc  OPA_decr_and_test_int
#define OPA_fetch_and_add_int_by_llsc  OPA_fetch_and_add_int
#define OPA_fetch_and_decr_int_by_llsc OPA_fetch_and_decr_int
#define OPA_fetch_and_incr_int_by_llsc OPA_fetch_and_incr_int
#define OPA_cas_ptr_by_llsc        OPA_cas_ptr
#define OPA_cas_int_by_llsc        OPA_cas_int
#define OPA_swap_ptr_by_llsc       OPA_swap_ptr
#define OPA_swap_int_by_llsc       OPA_swap_int


#include "opa_emulated.h"

#endif /* OPA_GCC_PPC_H_INCLUDED */
