/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2008 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef OPA_NT_INTRINSICS_H_INCLUDED
#define OPA_NT_INTRINSICS_H_INCLUDED

#define WIN32_LEAN_AND_MEAN
#include<windows.h>
#include<intrin.h>
#include "mpi.h"

/* OPA_int_t uses a long because the compiler intrinsics operate on
 * longs instead of ints. */
typedef struct { volatile long v;  } OPA_int_t;
typedef struct { void * volatile v; } OPA_ptr_t;

#define OPA_INT_T_INITIALIZER(val_) { (val_) }
#define OPA_PTR_T_INITIALIZER(val_) { (val_) }

#define OPA_write_barrier()      _WriteBarrier()
#define OPA_read_barrier()       _ReadBarrier()
#define OPA_read_write_barrier() _ReadWriteBarrier()
/* FIXME there mut be a more efficient way to implement this.  Is "asm {};"
 * sufficient? */
#define OPA_compiler_barrier()   _ReadWriteBarrier()

static _opa_inline int OPA_load_int(_opa_const OPA_int_t *ptr)
{
    return ((int)ptr->v);
}

static _opa_inline void OPA_store_int(OPA_int_t *ptr, int val)
{
    ptr->v = (long)val;
}

static _opa_inline void *OPA_load_ptr(_opa_const OPA_ptr_t *ptr)
{
    return ((void *)ptr->v);
}

static _opa_inline void OPA_store_ptr(OPA_ptr_t *ptr, void *val)
{
    ptr->v = val;
}

/* NOTE: these acquire/release operations have not been optimized, I just threw
 * down a full memory barrier.  Someone with more Windows expertise should feel
 * free to improve these (for Windows on x86/x86_64 these almost certainly only
 * need to be compiler barriers). */
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


static _opa_inline void OPA_add_int(OPA_int_t *ptr, int val)
{
    _InterlockedExchangeAdd(&(ptr->v), val);
}

static _opa_inline void OPA_incr_int(OPA_int_t *ptr)
{
    _InterlockedIncrement(&(ptr->v));
}

static _opa_inline void OPA_decr_int(OPA_int_t *ptr)
{
    _InterlockedDecrement(&(ptr->v));
}

static _opa_inline int OPA_decr_and_test_int(OPA_int_t *ptr)
{
    return (_InterlockedDecrement(&(ptr->v)) == 0);
}

static _opa_inline int OPA_fetch_and_add_int(OPA_int_t *ptr, int val)
{
    return ((int)_InterlockedExchangeAdd(&(ptr->v), (long)val));
}

static _opa_inline void *OPA_cas_ptr(OPA_ptr_t *ptr, void *oldv, void *newv)
{
#if (OPA_SIZEOF_VOID_P == 4)
    return ((void *)(LONG_PTR) _InterlockedCompareExchange((LONG volatile *)&(ptr->v),
                                                   (LONG)(LONG_PTR)newv,
                                                   (LONG)(LONG_PTR)oldv)
           );
#elif (OPA_SIZEOF_VOID_P == 8)
    return ((void *)(LONG_PTR)_InterlockedCompareExchange64((__int64 *)&(ptr->v),
                                                    (INT64)(LONG_PTR)newv,
                                                    (INT64)(LONG_PTR)oldv)
           );
#else
#error  "OPA_SIZEOF_VOID_P not valid"
#endif
}

static _opa_inline void *OPA_swap_ptr(OPA_ptr_t *ptr, void *val)
{
#if (OPA_SIZEOF_VOID_P == 4)
    return (void *)(LONG_PTR )_InterlockedExchange((LONG volatile *)&(ptr->v),
                                                    (LONG)(LONG_PTR)val);
#elif (OPA_SIZEOF_VOID_P == 8)
    return (void *)(LONG_PTR)_InterlockedExchange64((LONG64 volatile *)&(ptr->v),
                                                     (INT64)(LONG_PTR)val);
#else
#error  "OPA_SIZEOF_VOID_P not valid"
#endif
}

static _opa_inline int OPA_cas_int(OPA_int_t *ptr, int oldv, int newv)
{
    return _InterlockedCompareExchange((long *)&(ptr->v), newv, oldv);
}

static _opa_inline int OPA_swap_int(OPA_int_t *ptr, int val)
{
    return _InterlockedExchange(&(ptr->v), val);
}

/* Implement fetch_and_incr/decr using fetch_and_add (*_faa) */
#define OPA_fetch_and_incr_int_by_faa OPA_fetch_and_incr_int
#define OPA_fetch_and_decr_int_by_faa OPA_fetch_and_decr_int

#include "opa_emulated.h"

#endif /* defined(OPA_NT_INTRINSICS_H_INCLUDED) */
