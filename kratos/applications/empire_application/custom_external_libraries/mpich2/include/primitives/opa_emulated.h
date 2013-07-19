/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2008 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef OPA_EMULATED_H_INCLUDED
#define OPA_EMULATED_H_INCLUDED

/* Functions emulated using other atomics

   This header should be included at the bottom of any atomic
   primitives header that needs to implement an atomic op in terms of
   another atomic.
*/

/* Emulating using LL/SC */
#if defined(OPA_LL_SC_SUPPORTED)
static _opa_inline int OPA_fetch_and_add_int_by_llsc(OPA_int_t *ptr, int val)
{
    int prev;
    do {
        prev = OPA_LL_int(ptr);
    } while (!OPA_SC_int(ptr, prev + val));
    return prev;
}


static _opa_inline void OPA_add_int_by_llsc(OPA_int_t *ptr, int val)
{
    OPA_fetch_and_add_int_by_llsc(ptr, val);
}

static _opa_inline void OPA_incr_int_by_llsc(OPA_int_t *ptr)
{
    OPA_add_int_by_llsc(ptr, 1);
}

static _opa_inline void OPA_decr_int_by_llsc(OPA_int_t *ptr)
{
    OPA_add_int_by_llsc(ptr, -1);
}


static _opa_inline int OPA_fetch_and_decr_int_by_llsc(OPA_int_t *ptr)
{
    return OPA_fetch_and_add_int_by_llsc(ptr, -1);
}

static _opa_inline int OPA_fetch_and_incr_int_by_llsc(OPA_int_t *ptr)
{
    return OPA_fetch_and_add_int_by_llsc(ptr, 1);
}

static _opa_inline int OPA_decr_and_test_int_by_llsc(OPA_int_t *ptr)
{
    int prev = OPA_fetch_and_decr_int_by_llsc(ptr);
    return prev == 1;
}

static _opa_inline void *OPA_cas_ptr_by_llsc(OPA_ptr_t *ptr, void *oldv, void *newv)
{
    void *prev;
    do {
        prev = OPA_LL_ptr(ptr);
    } while (prev == oldv && !OPA_SC_ptr(ptr, newv));
    return prev;
}

static _opa_inline int OPA_cas_int_by_llsc(OPA_int_t *ptr, int oldv, int newv)
{
    int prev;
    do {
        prev = OPA_LL_int(ptr);
    } while (prev == oldv && !OPA_SC_int(ptr, newv));
    return prev;
}


static _opa_inline void *OPA_swap_ptr_by_llsc(OPA_ptr_t *ptr, void *val)
{
    void *prev;
    do {
        prev = OPA_LL_ptr(ptr);
    } while (!OPA_SC_ptr(ptr, val));
    return prev;
}

static _opa_inline int OPA_swap_int_by_llsc(OPA_int_t *ptr, int val)
{
    int prev;
    do {
        prev = OPA_LL_int(ptr);
    } while (!OPA_SC_int(ptr, val));
    return prev;
}

#endif /* OPA_LL_SC_SUPPORTED */

static _opa_inline int OPA_fetch_and_add_int_by_cas(OPA_int_t *ptr, int val)
{
    int cmp;
    int prev = OPA_load_int(ptr);

    do {
        cmp = prev;
        prev = OPA_cas_int(ptr, cmp, prev + val);
    } while (prev != cmp);

    return prev;
}

static _opa_inline int OPA_fetch_and_incr_int_by_faa(OPA_int_t *ptr)
{
    return OPA_fetch_and_add_int(ptr, 1);
}

static _opa_inline int OPA_fetch_and_decr_int_by_faa(OPA_int_t *ptr)
{
    return OPA_fetch_and_add_int(ptr, -1);
}

static _opa_inline int OPA_decr_and_test_int_by_fad(OPA_int_t *ptr)
{
    return OPA_fetch_and_decr_int(ptr) == 1;
}

static _opa_inline void OPA_add_int_by_faa(OPA_int_t *ptr, int val)
{
    OPA_fetch_and_add_int(ptr, val);
}

static _opa_inline int OPA_incr_int_by_faa(OPA_int_t *ptr)
{
    return OPA_fetch_and_add_int(ptr, 1);
}

static _opa_inline void OPA_incr_int_by_add(OPA_int_t *ptr)
{
    OPA_add_int(ptr, 1);
}

static _opa_inline void OPA_incr_int_by_fai(OPA_int_t *ptr)
{
    OPA_fetch_and_incr_int(ptr);
}

static _opa_inline int OPA_decr_int_by_faa(OPA_int_t *ptr)
{
    return OPA_fetch_and_add_int(ptr, -1);
}

static _opa_inline void OPA_decr_int_by_add(OPA_int_t *ptr)
{
    OPA_add_int(ptr, -1);
}

static _opa_inline void OPA_decr_int_by_fad(OPA_int_t *ptr)
{
    OPA_fetch_and_decr_int(ptr);
}


/* Swap using CAS */

static _opa_inline void *OPA_swap_ptr_by_cas(OPA_ptr_t *ptr, void *val)
{
    void *cmp;
    void *prev = OPA_load_ptr(ptr);

    do {
        cmp = prev;
        prev = OPA_cas_ptr(ptr, cmp, val);
    } while (prev != cmp);

    return prev;
}

static _opa_inline int OPA_swap_int_by_cas(OPA_int_t *ptr, int val)
{
    int cmp;
    int prev = OPA_load_int(ptr);

    do {
        cmp = prev;
        prev = OPA_cas_int(ptr, cmp, val);
    } while (prev != cmp);

    return prev;
}

#endif /* OPA_EMULATED_H_INCLUDED */
