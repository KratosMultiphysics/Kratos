/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2008 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef OPA_UTIL_H_INCLUDED
#define OPA_UTIL_H_INCLUDED

#define OPA_QUOTE(x_) OPA_QUOTE2(x_)
#define OPA_QUOTE2(x_) #x_

#if defined(OPA_HAVE_GCC_ATTRIBUTE)
#  define OPA_ATTRIBUTE(x_) __attribute__ (x_)
#else
#  define OPA_ATTRIBUTE(x_)
#endif

/* FIXME this just needs a total rework in general with an OPA_NDEBUG or similar. */
#define OPA_assert(expr_) do {} while (0)
#define OPA_assertp(expr_) do { if (!(expr_)) ++((int *)NULL) } while (0) /* SEGV intentionally */

/* A compile-time assertion macro.  It should cause a compilation error if (expr_) is false. */
#define OPA_COMPILE_TIME_ASSERT(expr_) \
    do { switch(0) { case 0: case (expr_): default: break; } } while (0)

#endif /* OPA_UTIL_H_INCLUDED */
