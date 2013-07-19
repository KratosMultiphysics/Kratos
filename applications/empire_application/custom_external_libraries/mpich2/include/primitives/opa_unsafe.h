/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2008 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef OPA_UNSAFE_H_INCLUDED
#define OPA_UNSAFE_H_INCLUDED

/* NOTE: These types intentionally do not use volatile in order to provide even
 * better performance.  If a good use case for having unsafe operations with
 * volatile types comes up, we can make two definitions that are conditional on
 * something sort of like NDEBUG. [goodell@ 2009-08-18] */
typedef struct { int v;  } OPA_int_t;
typedef struct { int *v; } OPA_ptr_t;

#define OPA_INT_T_INITIALIZER(val_) { (val_) }
#define OPA_PTR_T_INITIALIZER(val_) { (val_) }

/*
    Unsafe Primitives
    -----------------

    These are versions of the atomic primitives that emulate the single-threaded
    behavior of the primitives but do not attempt to provide any safety against
    concurrent use.  The primary use case for this implementation is to avoid
    overhead without code changes in client libraries and applications that have
    been compiled for single-threaded use only.  They might also be useful for
    testing (and testing of tests) in some cases.

    Thanks to Josh Haberman for inspiring this primitives implementation.
*/

static _opa_inline int OPA_load_int(const OPA_int_t *ptr)
{
    int retval;
    retval = ptr->v;
    return retval;
}

static _opa_inline void OPA_store_int(OPA_int_t *ptr, int val)
{
    ptr->v = val;
}

static _opa_inline void *OPA_load_ptr(const OPA_ptr_t *ptr)
{
    int *retval;
    retval = ptr->v;
    return retval;
}

static _opa_inline void OPA_store_ptr(OPA_ptr_t *ptr, void *val)
{
    ptr->v = val;
}

/* no barrier needed, not even a compiler barrier */
#define OPA_load_acquire_int(ptr_)       OPA_load_int((ptr_))
#define OPA_store_release_int(ptr_,val_) OPA_store_int((ptr_),(val_))
#define OPA_load_acquire_ptr(ptr_)       OPA_load_ptr((ptr_))
#define OPA_store_release_ptr(ptr_,val_) OPA_store_ptr((ptr_),(val_))

static _opa_inline void OPA_add_int(OPA_int_t *ptr, int val)
{
    ptr->v += val;
}

static _opa_inline void *OPA_cas_ptr(OPA_ptr_t *ptr, int *oldv, int *newv)
{
    int *prev;
    prev = ptr->v;
    if (prev == oldv) {
        ptr->v = newv;
    }
    return prev;
}

static _opa_inline int OPA_cas_int(OPA_int_t *ptr, int oldv, int newv)
{
    int prev;
    prev = ptr->v;
    if (prev == oldv) {
        ptr->v = newv;
    }
    return prev;
}

static _opa_inline int OPA_decr_and_test_int(OPA_int_t *ptr)
{
    int new_val;
    new_val = --(ptr->v);
    return (0 == new_val);
}

static _opa_inline void OPA_decr_int(OPA_int_t *ptr)
{
    --(ptr->v);
}

static _opa_inline int OPA_fetch_and_add_int(OPA_int_t *ptr, int val)
{
    int prev;
    prev = ptr->v;
    ptr->v += val;
    return prev;
}

static _opa_inline int OPA_fetch_and_decr_int(OPA_int_t *ptr)
{
    int prev;
    prev = ptr->v;
    --(ptr->v);
    return prev;
}

static _opa_inline int OPA_fetch_and_incr_int(OPA_int_t *ptr)
{
    int prev;
    prev = ptr->v;
    ++(ptr->v);
    return prev;
}

static _opa_inline void OPA_incr_int(OPA_int_t *ptr)
{
    ++(ptr->v);
}

static _opa_inline void *OPA_swap_ptr(OPA_ptr_t *ptr, void *val)
{
    int *prev;
    prev = ptr->v;
    ptr->v = val;
    return prev;
}

static _opa_inline int OPA_swap_int(OPA_int_t *ptr, int val)
{
    int prev;
    prev = ptr->v;
    ptr->v = val;
    return (int)prev;
}

/* null barriers */
#define OPA_write_barrier()      do {} while (0)
#define OPA_read_barrier()       do {} while (0)
#define OPA_read_write_barrier() do {} while (0)
#define OPA_compiler_barrier()   do {} while (0)

#endif /* !defined(OPA_UNSAFE_H_INCLUDED) */
