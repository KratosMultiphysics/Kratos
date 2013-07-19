/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*  
 *  (C) 2008 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#ifndef OPA_SUN_ATOMIC_OPS_H_INCLUDED
#define OPA_SUN_ATOMIC_OPS_H_INCLUDED

#include <atomic.h>

typedef struct { volatile uint_t v;  } OPA_int_t;
typedef struct { void * volatile v; } OPA_ptr_t;

#define OPA_INT_T_INITIALIZER(val_) { (val_) }
#define OPA_PTR_T_INITIALIZER(val_) { (val_) }

static _opa_inline int OPA_load_int(_opa_const OPA_int_t *ptr)
{
    return (int)ptr->v;
}

static _opa_inline void OPA_store_int(OPA_int_t *ptr, int val)
{
    ptr->v = (uint)val;
}

static _opa_inline void *OPA_load_ptr(_opa_const OPA_ptr_t *ptr)
{
    return ptr->v;
}

static _opa_inline void OPA_store_ptr(OPA_ptr_t *ptr, void *val)
{
    ptr->v = val;
}

static _opa_inline int OPA_load_acquire_int(_opa_const OPA_int_t *ptr)
{
    int tmp;
    tmp = ptr->v;
    membar_enter();
    return tmp;
}

static _opa_inline void OPA_store_release_int(OPA_int_t *ptr, int val)
{
    membar_exit();
    ptr->v = val;
}

static _opa_inline void *OPA_load_acquire_ptr(_opa_const OPA_ptr_t *ptr)
{
    void *tmp;
    tmp = ptr->v;
    membar_enter();
    return tmp;
}

static _opa_inline void OPA_store_release_ptr(OPA_ptr_t *ptr, void *val)
{
    membar_exit();
    ptr->v = val;
}


static _opa_inline void OPA_add_int(OPA_int_t *ptr, int val)
{
    atomic_add_int(&ptr->v, val);
}

static _opa_inline void OPA_incr_int(OPA_int_t *ptr)
{
    atomic_inc_uint(&ptr->v);
}

static _opa_inline void OPA_decr_int(OPA_int_t *ptr)
{
    atomic_dec_uint(&ptr->v);
}


static _opa_inline int OPA_decr_and_test_int(OPA_int_t *ptr)
{
    return atomic_dec_uint_nv(&ptr->v) == 0;    
}


static _opa_inline int OPA_fetch_and_add_int(OPA_int_t *ptr, int val)
{
    return (int)atomic_add_int_nv(&ptr->v, val) - val;
}

static _opa_inline int OPA_fetch_and_decr_int(OPA_int_t *ptr)
{
    return (int)atomic_dec_uint_nv(&ptr->v) + 1;
}

static _opa_inline int OPA_fetch_and_incr_int(OPA_int_t *ptr)
{
    return (int)atomic_inc_uint_nv(&ptr->v) - 1;
}


static _opa_inline void *OPA_cas_ptr(OPA_ptr_t *ptr, void *oldv, void *newv)
{
    return atomic_cas_ptr(ptr, oldv, newv);
}

static _opa_inline int OPA_cas_int(OPA_int_t *ptr, int oldv, int newv)
{
    return (int)atomic_cas_uint(&ptr->v, (uint_t)oldv, (uint_t)newv);
}


static _opa_inline void *OPA_swap_ptr(OPA_ptr_t *ptr, void *val)
{
    return atomic_swap_ptr(ptr, val);
}

static _opa_inline int OPA_swap_int(OPA_int_t *ptr, int val)
{
    return (int)atomic_swap_uint(&ptr->v, (uint_t) val);
}


#define OPA_write_barrier()      membar_producer()
#define OPA_read_barrier()       membar_consumer()
#define OPA_read_write_barrier() do { membar_consumer(); membar_producer(); } while (0)

/* is this portable enough? */
#define OPA_compiler_barrier()   __asm__ __volatile__  ( ""  ::: "memory" )

#endif /* OPA_SUN_ATOMIC_OPS_H_INCLUDED */
