/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*  
 *  (C) 2008 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef OPA_PRIMITIVES_H_INCLUDED
#define OPA_PRIMITIVES_H_INCLUDED

#include "opa_config.h"
#include "opa_util.h"

/* Clean up some of the opa_config.h definitions.  This is a consequence
   of using the AX_PREFIX_CONFIG_H macro. Autoconf won't define inline
   or _opa_inline when a real "inline" works.  Since we are
   unconditionally using _opa_inline we must define it ourselves in this
   case. */
#ifndef _opa_inline
#define _opa_inline inline
#endif
#ifndef _opa_restrict
#define _opa_restrict restrict
#endif
#ifndef _opa_const
#define _opa_const const
#endif

/*
   Primitive atomic functions
   --------------------------

   The included file is responsible for defining the types of OPA_int_t and
   OPA_ptr_t as well as a set of functions for operating on these
   types.  If you have the following declaration:

       OPA_int_t atomic_var;

   Then in order for the emulation functions to compile, the underlying value of
   atomic_var should be accessible via:

       atomic_var.v;

   The same goes for OPA_ptr_t.

   The atomic functions that must be ported for each architecture: 

   static _opa_inline int   OPA_load_int(_opa_const OPA_int_t *ptr);
   static _opa_inline void  OPA_store_int(OPA_int_t *ptr, int val);
   static _opa_inline void *OPA_load_ptr(_opa_const OPA_ptr_t *ptr);
   static _opa_inline void  OPA_store_ptr(OPA_ptr_t *ptr, void *val);

   static _opa_inline void OPA_add_int(OPA_int_t *ptr, int val);
   static _opa_inline void OPA_incr_int(OPA_int_t *ptr);
   static _opa_inline void OPA_decr_int(OPA_int_t *ptr);

   static _opa_inline int OPA_decr_and_test_int(OPA_int_t *ptr);
   static _opa_inline int OPA_fetch_and_add_int(OPA_int_t *ptr, int val);
   static _opa_inline int OPA_fetch_and_incr_int(OPA_int_t *ptr);
   static _opa_inline int OPA_fetch_and_decr_int(OPA_int_t *ptr);

   static _opa_inline void *OPA_cas_ptr(OPA_ptr_t *ptr, void *oldv, void *newv);
   static _opa_inline int   OPA_cas_int(OPA_int_t *ptr, int oldv, int newv);

   static _opa_inline void *OPA_swap_ptr(OPA_ptr_t *ptr, void *val);
   static _opa_inline int   OPA_swap_int(OPA_int_t *ptr, int val);

   // (the memory barriers may be macros instead of inline functions)
   static _opa_inline void OPA_write_barrier();
   static _opa_inline void OPA_read_barrier();
   static _opa_inline void OPA_read_write_barrier();

   // Loads and stores with memory ordering guarantees (also may be macros):
   static _opa_inline int   OPA_load_acquire_int(_opa_const OPA_int_t *ptr);
   static _opa_inline void  OPA_store_release_int(OPA_int_t *ptr, int val);
   static _opa_inline void *OPA_load_acquire_ptr(_opa_const OPA_ptr_t *ptr);
   static _opa_inline void  OPA_store_release_ptr(OPA_ptr_t *ptr, void *val);

   // Compiler barrier, only preventing compiler reordering, *not* CPU
   // reordering (may be a macro):
   static _opa_inline void OPA_compiler_barrier();

   // The following need to be ported only for architectures supporting LL/SC:
   static _opa_inline int OPA_LL_int(OPA_int_t *ptr);
   static _opa_inline int OPA_SC_int(OPA_int_t *ptr, int val);
   static _opa_inline void *OPA_LL_ptr(OPA_ptr_t *ptr);
   static _opa_inline int OPA_SC_ptr(OPA_ptr_t *ptr, void *val);

   // Additionally, the following initializer macros must be defined:
   #define OPA_INT_T_INITIALIZER(val_) ...
   #define OPA_PTR_T_INITIALIZER(val_) ...

   // They should be useable as C89 static initializers like so:

   struct { int x; OPA_int_t y; OPA_ptr_t z; } foo = { 35, OPA_INT_T_INITIALIZER(1), OPA_PTR_T_INITIALIZER(NULL) };
*/

/* Include the appropriate header for the architecture */
#if defined(OPA_USE_UNSAFE_PRIMITIVES)
/* comes first to permit user overrides in the style of NDEBUG */
#include "primitives/opa_unsafe.h"
#elif   defined(OPA_HAVE_GCC_AND_POWERPC_ASM)
#include "primitives/opa_gcc_ppc.h"
#elif defined(OPA_HAVE_GCC_X86_32_64)
#include "primitives/opa_gcc_intel_32_64.h"
#elif defined(OPA_HAVE_GCC_X86_32_64_P3)
#include "primitives/opa_gcc_intel_32_64_p3.h"
#elif defined(OPA_HAVE_GCC_AND_IA64_ASM)
#include "primitives/opa_gcc_ia64.h"
#elif defined(OPA_HAVE_GCC_AND_SICORTEX_ASM)
#include "primitives/opa_gcc_sicortex.h"
#elif defined(OPA_HAVE_GCC_INTRINSIC_ATOMICS)
#include "primitives/opa_gcc_intrinsics.h"
#elif defined(OPA_HAVE_SUN_ATOMIC_OPS)
#include "primitives/opa_sun_atomic_ops.h"
#elif defined(OPA_HAVE_NT_INTRINSICS)
#include "primitives/opa_nt_intrinsics.h"
#elif defined(OPA_USE_LOCK_BASED_PRIMITIVES)
#include "primitives/opa_by_lock.h"
#else
#error no primitives implementation specified
#endif

/*
    This routine is needed because the MPIU_THREAD_XXX_CS_{ENTER,EXIT} macros do
    not provide synchronization across multiple processes, only across multiple
    threads within a process.  In order to safely emulate atomic operations on a
    shared memory region, we need a shared memory backed lock mechanism.

    This routine must be called by any subsystem that intends to use the atomic
    abstractions if the cpp directive OPA_USE_LOCK_BASED_PRIMITIVES is defined.  It must
    be called exactly once by _all_ processes, not just a single leader.  This
    function will initialize the contents of the lock variable if the caller
    specifies (isLeader==true).  Note that multiple initialization is forbidden
    by several lock implementations, especially pthreads.

    Inputs:
      shm_lock - A pointer to an allocated piece of shared memory that can hold
                 a mutex (e.g., pthread_mutex_t).  This is not portable to
                 non-pthreads systems at this time.
      isLeader - This boolean value should be set to true for exactly one
                 thread/process of the group that calls this function.
*/
#if defined(OPA_HAVE_PTHREAD_H)
#  include <pthread.h>
typedef pthread_mutex_t OPA_emulation_ipl_t;
int OPA_Interprocess_lock_init(OPA_emulation_ipl_t *shm_lock, int isLeader);
#endif


/* FIXME This should probably be pushed down into the platform-specific headers. */
#if defined(OPA_HAVE_SCHED_YIELD)
#  include <sched.h>
#  define OPA_busy_wait() sched_yield()
#else
#  define OPA_busy_wait() do { } while (0)
#endif

#endif /* defined(OPA_PRIMITIVES_H_INCLUDED) */
