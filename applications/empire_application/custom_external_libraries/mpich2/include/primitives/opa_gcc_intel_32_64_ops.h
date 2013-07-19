/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2008 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef OPA_GCC_INTEL_32_64_OPS_H_INCLUDED
#define OPA_GCC_INTEL_32_64_OPS_H_INCLUDED

#ifndef OPA_SIZEOF_INT
#error OPA_SIZEOF_INT is not defined
#endif

/* Set OPA_SS (Size Suffix) which is to be appended to asm ops for
   specifying 4 or 8 byte operands */
#if OPA_SIZEOF_INT == 4
#define OPA_SS "l"
#elif OPA_SIZEOF_INT == 8
#define OPA_SS "q"
#else
#error OPA_SIZEOF_INT is not 4 or 8
#endif

/* XXX DJG FIXME do we need to align these? */
typedef struct { volatile int v;    } OPA_int_t;
typedef struct { void * volatile v; } OPA_ptr_t;

#define OPA_INT_T_INITIALIZER(val_) { (val_) }
#define OPA_PTR_T_INITIALIZER(val_) { (val_) }

/* Aligned loads and stores are atomic on x86(-64). */
static _opa_inline int OPA_load_int(_opa_const OPA_int_t *ptr)
{
    return ptr->v;
}

/* Aligned loads and stores are atomic on x86(-64). */
static _opa_inline void OPA_store_int(OPA_int_t *ptr, int val)
{
    ptr->v = val;
}

/* Aligned loads and stores are atomic on x86(-64). */
static _opa_inline void *OPA_load_ptr(_opa_const OPA_ptr_t *ptr)
{
    return ptr->v;
}

/* Aligned loads and stores are atomic on x86(-64). */
static _opa_inline void OPA_store_ptr(OPA_ptr_t *ptr, void *val)
{
    ptr->v = val;
}

/* NOTE-X86-1 x86 and x86_64 processors execute instructions such that their
 * effects are visible in issue order, with one exception:  earlier stores may
 * pass later loads.  Also, non-temporal stores (with the "movnt*" instructions)
 * do not follow these ordering constraints.  So all normal load/store
 * instructions on x86(_64) already have acquire/release semantics.  We just
 * have to make sure that the compiler does not reorder the instructions. */

static _opa_inline int OPA_load_acquire_int(_opa_const OPA_int_t *ptr)
{
    int tmp;
    tmp = ptr->v;
    OPA_compiler_barrier(); /* NOTE-X86-1 */
    return tmp;
}

static _opa_inline void OPA_store_release_int(OPA_int_t *ptr, int val)
{
    OPA_compiler_barrier(); /* NOTE-X86-1 */
    ptr->v = val;
}

static _opa_inline void *OPA_load_acquire_ptr(_opa_const OPA_ptr_t *ptr)
{
    void *tmp;
    tmp = ptr->v;
    OPA_compiler_barrier(); /* NOTE-X86-1 */
    return tmp;
}

static _opa_inline void OPA_store_release_ptr(OPA_ptr_t *ptr, void *val)
{
    OPA_compiler_barrier(); /* NOTE-X86-1 */
    ptr->v = val;
}


static _opa_inline void OPA_add_int(OPA_int_t *ptr, int val)
{
    __asm__ __volatile__ ("lock ; add"OPA_SS" %1,%0"
                          :"=m" (ptr->v)
                          :"ir" (val), "m" (ptr->v));
    return;
}

static _opa_inline void OPA_incr_int(OPA_int_t *ptr)
{
    __asm__ __volatile__ ("lock ; inc"OPA_SS" %0"
                          :"=m" (ptr->v)
                          :"m" (ptr->v));
    return;
}

static _opa_inline void OPA_decr_int(OPA_int_t *ptr)
{
    __asm__ __volatile__ ("lock ; dec"OPA_SS" %0"
                          :"=m" (ptr->v)
                          :"m" (ptr->v));
    return;
}


static _opa_inline int OPA_decr_and_test_int(OPA_int_t *ptr)
{
    char result;
    __asm__ __volatile__ ("lock ; dec"OPA_SS" %0; setz %1"
                          :"=m" (ptr->v), "=q" (result)
                          :"m" (ptr->v));
    return result;
}

static _opa_inline int OPA_fetch_and_add_int(OPA_int_t *ptr, int val)
{
    __asm__ __volatile__ ("lock ; xadd %0,%1"
                          : "=r" (val), "=m" (ptr->v)
                          :  "0" (val),  "m" (ptr->v));
    return val;
}

#define OPA_fetch_and_incr_int_by_faa OPA_fetch_and_incr_int
#define OPA_fetch_and_decr_int_by_faa OPA_fetch_and_decr_int


static _opa_inline void *OPA_cas_ptr(OPA_ptr_t *ptr, void *oldv, void *newv)
{
    void *prev;
    __asm__ __volatile__ ("lock ; cmpxchg %3,%4"
                          : "=a" (prev), "=m" (ptr->v)
                          : "0" (oldv), "q" (newv), "m" (ptr->v));
    return prev;
}

static _opa_inline int OPA_cas_int(OPA_int_t *ptr, int oldv, int newv)
{
    int prev;
    __asm__ __volatile__ ("lock ; cmpxchg %3,%4"
                          : "=a" (prev), "=m" (ptr->v)
                          : "0" (oldv), "q" (newv), "m" (ptr->v));
    return prev;
}

static _opa_inline void *OPA_swap_ptr(OPA_ptr_t *ptr, void *val)
{
    __asm__ __volatile__ ("xchg %0,%1"
                          :"=r" (val), "=m" (ptr->v)
                          : "0" (val),  "m" (ptr->v));
    return val;
}

static _opa_inline int OPA_swap_int(OPA_int_t *ptr, int val)
{
    __asm__ __volatile__ ("xchg %0,%1"
                          :"=r" (val), "=m" (ptr->v)
                          : "0" (val),  "m" (ptr->v));
    return val;
}

#undef OPA_SS

#endif /* OPA_GCC_INTEL_32_64_OPS_H_INCLUDED */
