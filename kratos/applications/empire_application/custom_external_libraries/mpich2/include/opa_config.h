#ifndef _SRC_OPA_CONFIG_H
#define _SRC_OPA_CONFIG_H 1
 
/* src/opa_config.h. Generated automatically at end of configure. */
/* src/config.h.  Generated from config.h.in by configure.  */
/* src/config.h.in.  Generated from configure.ac by autoheader.  */

/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2008 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


/* define if lock-based emulation was explicitly requested at configure time
   via --with-atomic-primitives=no */
/* #undef EXPLICIT_EMULATION */

/* Define to 1 if you have the <atomic.h> header file. */
/* #undef HAVE_ATOMIC_H */

/* Define to 1 if you have the <dlfcn.h> header file. */
#ifndef OPA_HAVE_DLFCN_H 
#define OPA_HAVE_DLFCN_H  1 
#endif

/* define to 1 if we have support for gcc ia64 primitives */
/* #undef HAVE_GCC_AND_IA64_ASM */

/* define to 1 if we have support for gcc PowerPC atomics */
/* #undef HAVE_GCC_AND_POWERPC_ASM */

/* define to 1 if we have support for gcc SiCortex atomics */
/* #undef HAVE_GCC_AND_SICORTEX_ASM */

/* Define if GNU __attribute__ is supported */
#ifndef OPA_HAVE_GCC_ATTRIBUTE 
#define OPA_HAVE_GCC_ATTRIBUTE  1 
#endif

/* define to 1 if we have support for gcc atomic intrinsics */
#ifndef OPA_HAVE_GCC_INTRINSIC_ATOMICS 
#define OPA_HAVE_GCC_INTRINSIC_ATOMICS  1 
#endif

/* define to 1 if we have support for gcc x86/x86_64 primitives */
#ifndef OPA_HAVE_GCC_X86_32_64 
#define OPA_HAVE_GCC_X86_32_64  1 
#endif

/* define to 1 if we have support for gcc x86 primitives for pre-Pentium 4 */
#ifndef OPA_HAVE_GCC_X86_32_64_P3 
#define OPA_HAVE_GCC_X86_32_64_P3  1 
#endif

/* Define to 1 if you have the <intrin.h> header file. */
/* #undef HAVE_INTRIN_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#ifndef OPA_HAVE_INTTYPES_H 
#define OPA_HAVE_INTTYPES_H  1 
#endif

/* Define to 1 if you have the `pthread' library (-lpthread). */
#ifndef OPA_HAVE_LIBPTHREAD 
#define OPA_HAVE_LIBPTHREAD  1 
#endif

/* Define to 1 if you have the <memory.h> header file. */
#ifndef OPA_HAVE_MEMORY_H 
#define OPA_HAVE_MEMORY_H  1 
#endif

/* define to 1 if we have support for Windows NT atomic intrinsics */
/* #undef HAVE_NT_INTRINSICS */

/* Define to 1 if you have the <pthread.h> header file. */
#ifndef OPA_HAVE_PTHREAD_H 
#define OPA_HAVE_PTHREAD_H  1 
#endif

/* Define to 1 if you have the `pthread_yield' function. */
#ifndef OPA_HAVE_PTHREAD_YIELD 
#define OPA_HAVE_PTHREAD_YIELD  1 
#endif

/* Define to 1 if you have the `sched_yield' function. */
/* #undef HAVE_SCHED_YIELD */

/* Define to 1 if you have the <stddef.h> header file. */
#ifndef OPA_HAVE_STDDEF_H 
#define OPA_HAVE_STDDEF_H  1 
#endif

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef OPA_HAVE_STDINT_H 
#define OPA_HAVE_STDINT_H  1 
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef OPA_HAVE_STDLIB_H 
#define OPA_HAVE_STDLIB_H  1 
#endif

/* Define if strict checking of atomic operation fairness is desired */
/* #undef HAVE_STRICT_FAIRNESS_CHECKS */

/* Define to 1 if you have the <strings.h> header file. */
#ifndef OPA_HAVE_STRINGS_H 
#define OPA_HAVE_STRINGS_H  1 
#endif

/* Define to 1 if you have the <string.h> header file. */
#ifndef OPA_HAVE_STRING_H 
#define OPA_HAVE_STRING_H  1 
#endif

/* define to 1 if we have support for Sun atomic operations library */
/* #undef HAVE_SUN_ATOMIC_OPS */

/* Define to 1 if you have the <sys/stat.h> header file. */
#ifndef OPA_HAVE_SYS_STAT_H 
#define OPA_HAVE_SYS_STAT_H  1 
#endif

/* Define to 1 if you have the <sys/types.h> header file. */
#ifndef OPA_HAVE_SYS_TYPES_H 
#define OPA_HAVE_SYS_TYPES_H  1 
#endif

/* Define to 1 if you have the <unistd.h> header file. */
#ifndef OPA_HAVE_UNISTD_H 
#define OPA_HAVE_UNISTD_H  1 
#endif

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#ifndef OPA_LT_OBJDIR 
#define OPA_LT_OBJDIR  ".libs/" 
#endif

/* define to the maximum number of simultaneous threads */
#ifndef OPA_MAX_NTHREADS 
#define OPA_MAX_NTHREADS  100 
#endif

/* Define to 1 if assertions should be disabled. */
/* #undef NDEBUG */

/* Name of package */
#ifndef OPA_PACKAGE 
#define OPA_PACKAGE  "openpa" 
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef OPA_PACKAGE_BUGREPORT 
#define OPA_PACKAGE_BUGREPORT  "https://trac.mcs.anl.gov/projects/openpa/newticket" 
#endif

/* Define to the full name of this package. */
#ifndef OPA_PACKAGE_NAME 
#define OPA_PACKAGE_NAME  "OpenPA" 
#endif

/* Define to the full name and version of this package. */
#ifndef OPA_PACKAGE_STRING 
#define OPA_PACKAGE_STRING  "OpenPA 1.0.3" 
#endif

/* Define to the one symbol short name of this package. */
#ifndef OPA_PACKAGE_TARNAME 
#define OPA_PACKAGE_TARNAME  "openpa" 
#endif

/* Define to the home page for this package. */
#ifndef OPA_PACKAGE_URL 
#define OPA_PACKAGE_URL  "" 
#endif

/* Define to the version of this package. */
#ifndef OPA_PACKAGE_VERSION 
#define OPA_PACKAGE_VERSION  "1.0.3" 
#endif

/* The size of `int', as computed by sizeof. */
#ifndef OPA_SIZEOF_INT 
#define OPA_SIZEOF_INT  4 
#endif

/* The size of `void *', as computed by sizeof. */
#ifndef OPA_SIZEOF_VOID_P 
#define OPA_SIZEOF_VOID_P  8 
#endif

/* Define to 1 if you have the ANSI C header files. */
#ifndef OPA_STDC_HEADERS 
#define OPA_STDC_HEADERS  1 
#endif

/* define to 1 to force using lock-based atomic primitives */
/* #undef USE_LOCK_BASED_PRIMITIVES */

/* define to 1 if unsafe (non-atomic) primitives should be used */
/* #undef USE_UNSAFE_PRIMITIVES */

/* Version number of package */
#ifndef OPA_VERSION 
#define OPA_VERSION  "1.0.3" 
#endif

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported directly.  */
#ifndef _opa_restrict 
#define _opa_restrict  __restrict 
#endif
/* Work around a bug in Sun C++: it does not support _Restrict or
   __restrict__, even though the corresponding Sun C compiler ends up with
   "#define restrict _Restrict" or "#define restrict __restrict__" in the
   previous line.  Perhaps some future version of Sun C++ will work with
   restrict; if so, hopefully it defines __RESTRICT like Sun C does.  */
#if defined __SUNPRO_CC && !defined __RESTRICT
# define _Restrict
# define __restrict__
#endif


 
/* once: _SRC_OPA_CONFIG_H */
#endif
