/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2008 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef OPA_GCC_INTEL_32_64_P3BARRIER_H_INCLUDED
#define OPA_GCC_INTEL_32_64_P3BARRIER_H_INCLUDED

#define OPA_compiler_barrier() __asm__ __volatile__ ( "" ::: "memory" )

/* For all regular memory (write-back cacheable, not driver/graphics
 * memory), there is only one general ordering relaxation permitted by
 * x86/x86_64 processors: earlier stores may be retired after later
 * stores.  The "clflush" and "movnt*" instructions also don't follow
 * general ordering constraints, although any code using these
 * instructions should be responsible for ensuring proper ordering
 * itself.  So our read and write barriers may be implemented as simple
 * compiler barriers. */
#define OPA_write_barrier() OPA_compiler_barrier()
#define OPA_read_barrier()  OPA_compiler_barrier()

/* Some Pentium III and earlier processors have store barriers
   (sfence), but no full barrier or load barriers.  Some other very
   recent x86-like processors don't seem to have sfence either.  A
   lock-prefixed atomic instruction (operating on any memory) behaves as
   a full memory barrier, so we use that here. */
static inline void OPA_read_write_barrier(void)
{
    int a = 0;
    __asm__ __volatile__ ("lock orl $0, %0" : "+m" (a) : /*no outputs*/ : "memory");
}

#endif /* OPA_GCC_INTEL_32_64_P3BARRIER_H_INCLUDED */
