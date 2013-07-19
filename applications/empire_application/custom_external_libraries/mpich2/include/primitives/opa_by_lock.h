/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2008 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef OPA_BY_LOCK_H_INCLUDED
#define OPA_BY_LOCK_H_INCLUDED

/* FIXME For now we rely on pthreads for our IPC locks.  This is fairly
   portable, although it is obviously not 100% portable.  Some day when we
   refactor the OPA_Process_locks code we should be able to use that again. */
#if defined(OPA_HAVE_PTHREAD_H)
#include <pthread.h>

/* defined in opa_primitives.c */
extern pthread_mutex_t *OPA_emulation_lock;

/* FIXME these make less sense now that OPA is not inside of MPICH2.  Is there a
   simpler name/scheme that could be used here instead? [goodell@ 2009-02-19] */
#define OPA_IPC_SINGLE_CS_ENTER(msg)          \
    do {                                        \
        OPA_assert(OPA_emulation_lock);    \
        pthread_mutex_lock(OPA_emulation_lock);     \
    } while (0)

#define OPA_IPC_SINGLE_CS_EXIT(msg)           \
    do {                                        \
        OPA_assert(OPA_emulation_lock);    \
        pthread_mutex_unlock(OPA_emulation_lock);   \
    } while (0)

typedef struct { volatile int v;  } OPA_int_t;
typedef struct { int * volatile v; } OPA_ptr_t;

#define OPA_INT_T_INITIALIZER(val_) { (val_) }
#define OPA_PTR_T_INITIALIZER(val_) { (val_) }

/*
    Emulated atomic primitives
    --------------------------

    These are versions of the atomic primitives that emulate the proper behavior
    via the use of an inter-process lock.  For more information on their
    individual behavior, please see the comment on the corresponding top level
    function.

    In general, these emulated primitives should _not_ be used.  Most algorithms
    can be more efficiently implemented by putting most or all of the algorithm
    inside of a single critical section.  These emulated primitives exist to
    ensure that there is always a fallback if no machine-dependent version of a
    particular operation has been defined.  They also serve as a very readable
    reference for the exact semantics of our OPA_* ops.
*/

static _opa_inline int OPA_load_int(_opa_const OPA_int_t *ptr)
{
    int retval;
    OPA_IPC_SINGLE_CS_ENTER("atomic_add");
    retval = ptr->v;
    OPA_IPC_SINGLE_CS_EXIT("atomic_add");
    return retval;
}

static _opa_inline void OPA_store_int(OPA_int_t *ptr, int val)
{
    OPA_IPC_SINGLE_CS_ENTER("atomic_add");
    ptr->v = val;
    OPA_IPC_SINGLE_CS_EXIT("atomic_add");
}

static _opa_inline void *OPA_load_ptr(_opa_const OPA_ptr_t *ptr)
{
    int *retval;
    OPA_IPC_SINGLE_CS_ENTER("atomic_add");
    retval = ptr->v;
    OPA_IPC_SINGLE_CS_EXIT("atomic_add");
    return retval;
}

static _opa_inline void OPA_store_ptr(OPA_ptr_t *ptr, void *val)
{
    OPA_IPC_SINGLE_CS_ENTER("atomic_add");
    ptr->v = val;
    OPA_IPC_SINGLE_CS_EXIT("atomic_add");
}

/* normal loads/stores are fully ordered, so just use them */
#define OPA_load_acquire_int(ptr_)       OPA_load_int((ptr_))
#define OPA_store_release_int(ptr_,val_) OPA_store_int((ptr_),(val_))
#define OPA_load_acquire_ptr(ptr_)       OPA_load_ptr((ptr_))
#define OPA_store_release_ptr(ptr_,val_) OPA_store_ptr((ptr_),(val_))

static _opa_inline void OPA_add_int(OPA_int_t *ptr, int val)
{
    OPA_IPC_SINGLE_CS_ENTER("atomic_add");
    ptr->v += val;
    OPA_IPC_SINGLE_CS_EXIT("atomic_add");
}

static _opa_inline void *OPA_cas_ptr(OPA_ptr_t *ptr, int *oldv, int *newv)
{
    int *prev;
    OPA_IPC_SINGLE_CS_ENTER("atomic_cas");
    prev = ptr->v;
    if (prev == oldv) {
        ptr->v = newv;
    }
    OPA_IPC_SINGLE_CS_EXIT("atomic_cas");
    return prev;
}

static _opa_inline int OPA_cas_int(OPA_int_t *ptr, int oldv, int newv)
{
    int prev;
    OPA_IPC_SINGLE_CS_ENTER("atomic_cas");
    prev = ptr->v;
    if (prev == oldv) {
        ptr->v = newv;
    }
    OPA_IPC_SINGLE_CS_EXIT("atomic_cas");
    return prev;
}

static _opa_inline int OPA_decr_and_test_int(OPA_int_t *ptr)
{
    int new_val;
    OPA_IPC_SINGLE_CS_ENTER("atomic_decr_and_test");
    new_val = --(ptr->v);
    OPA_IPC_SINGLE_CS_EXIT("atomic_decr_and_test");
    return (0 == new_val);
}

static _opa_inline void OPA_decr_int(OPA_int_t *ptr)
{
    OPA_IPC_SINGLE_CS_ENTER("atomic_decr");
    --(ptr->v);
    OPA_IPC_SINGLE_CS_EXIT("atomic_decr");
}

static _opa_inline int OPA_fetch_and_add_int(OPA_int_t *ptr, int val)
{
    int prev;
    OPA_IPC_SINGLE_CS_ENTER("atomic_fetch_and_add");
    prev = ptr->v;
    ptr->v += val;
    OPA_IPC_SINGLE_CS_EXIT("atomic_fetch_and_add");
    return prev;
}

static _opa_inline int OPA_fetch_and_decr_int(OPA_int_t *ptr)
{
    int prev;
    OPA_IPC_SINGLE_CS_ENTER("atomic_fetch_and_decr");
    prev = ptr->v;
    --(ptr->v);
    OPA_IPC_SINGLE_CS_EXIT("atomic_fetch_and_decr");
    return prev;
}

static _opa_inline int OPA_fetch_and_incr_int(OPA_int_t *ptr)
{
    int prev;
    OPA_IPC_SINGLE_CS_ENTER("atomic_fetch_and_incr");
    prev = ptr->v;
    ++(ptr->v);
    OPA_IPC_SINGLE_CS_EXIT("atomic_fetch_and_incr");
    return prev;
}

static _opa_inline void OPA_incr_int(OPA_int_t *ptr)
{
    OPA_IPC_SINGLE_CS_ENTER("atomic_incr");
    ++(ptr->v);
    OPA_IPC_SINGLE_CS_EXIT("atomic_incr");
}

static _opa_inline void *OPA_swap_ptr(OPA_ptr_t *ptr, void *val)
{
    int *prev;
    OPA_IPC_SINGLE_CS_ENTER("atomic_swap_ptr");
    prev = ptr->v;
    ptr->v = val;
    OPA_IPC_SINGLE_CS_EXIT("atomic_swap_ptr");
    return prev;
}

static _opa_inline int OPA_swap_int(OPA_int_t *ptr, int val)
{
    int prev;
    OPA_IPC_SINGLE_CS_ENTER("atomic_swap_int");
    prev = ptr->v;
    ptr->v = val;
    OPA_IPC_SINGLE_CS_EXIT("atomic_swap_int");
    return (int)prev;
}

/* lock/unlock provides barrier */
#define OPA_write_barrier()      do {} while (0)
#define OPA_read_barrier()       do {} while (0)
#define OPA_read_write_barrier() do {} while (0)


#endif /* defined(OPA_HAVE_PTHREAD_H) */
#endif /* !defined(OPA_BY_LOCK_H_INCLUDED) */
