/*
 * Copyright (c) 2009-2010 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 */

#ifndef __PPC
#warning This file is only for PowerPC
#endif  /* __ PPC */

static inline void pastix_mfence( void )
{
    __asm__ __volatile__ ("lwsync\n\t":::"memory");
}

static inline int pastix_atomic_bor_32b( volatile uint32_t* location,
                                          uint32_t mask )
{
#if !defined(__IBMC__)
   int32_t old, t;

   __asm__ __volatile__(
                        "1:   lwarx   %0,  0, %3   \n\t"
                        "     or      %1, %0, %2   \n\t"
                        "     stwcx.  %1,  0, %3   \n\t"
                        "     bne-    1b           \n\t"
                        : "=&r" (old), "=&r" (t)
                        : "r" (mask), "r" (location)
                        : "cc", "memory");

   return t;
#else
   return mask | __fetch_and_or(location, mask);
#endif  /* !defined(__IBMC__) */
}

static inline int pastix_atomic_band_32b( volatile uint32_t* location,
                                          uint32_t mask )
{
#if !defined(__IBMC__)
   int32_t old, t;

   __asm__ __volatile__(
                        "1:   lwarx   %0,  0, %3   \n\t"
                        "     andc    %1, %0, %2   \n\t"
                        "     stwcx.  %1,  0, %3   \n\t"
                        "     bne-    1b           \n\t"
                        : "=&r" (old), "=&r" (t)
                        : "r" (mask), "r" (location)
                        : "cc", "memory");

   return t;
#else
   return mask & __fetch_and_and(location, mask);
#endif  /* !defined(__IBMC__) */
}

static inline int pastix_atomic_cas_32b( volatile uint32_t* location,
                                          uint32_t old_value,
                                          uint32_t new_value )
{
#if !defined(__IBMC__)
   int32_t ret;

   __asm__ __volatile__ (
                         "1: lwarx   %0, 0, %2  \n\t"
                         "   cmpw    0, %0, %3  \n\t"
                         "   bne-    2f         \n\t"
                         "   stwcx.  %4, 0, %2  \n\t"
                         "   bne-    1b         \n\t"
                         "2:"
                         : "=&r" (ret), "=m" (*location)
                         : "r" (location), "r" (old_value), "r" (new_value), "m" (*location)
                         : "cr0", "memory");

   return (ret == old_value);
#else
   return __compare_and_swap((volatile int*)location, (int*)&old_value, (int)new_value);
#endif  /* !defined(__IBMC__) */
}

static inline int pastix_atomic_cas_64b( volatile uint64_t* location,
                                          uint64_t old_value,
                                          uint64_t new_value )
{
#if !defined(__IBMC__)
   int64_t ret;

   __asm__ __volatile__ (
                         "1: ldarx   %0, 0, %2  \n\t"
                         "   cmpd    0, %0, %3  \n\t"
                         "   bne-    2f         \n\t"
                         "   stdcx.  %4, 0, %2  \n\t"
                         "   bne-    1b         \n\t"
                         "2:"
                         : "=&r" (ret), "=m" (*location)
                         : "r" (location), "r" (old_value), "r" (new_value), "m" (*location)
                         : "cr0", "memory");

   return (ret == old_value);
#else
   return __compare_and_swaplp((volatile long*)location, (long*)&old_value, (long)new_value);
#endif  /* !defined(__IBMC__) */
}

#define PASTIX_ATOMIC_HAS_ATOMIC_INC_32B
static inline int32_t pastix_atomic_inc_32b( volatile int32_t *location )
{
#if !defined(__IBMC__)
   int32_t t;

   __asm__ __volatile__(
                        "1:   lwarx   %0, 0, %1    \n\t"
                        "     addic   %0, %0, 1    \n\t"
                        "     stwcx.  %0, 0, %1    \n\t"
                        "     bne-    1b           \n\t"
                        : "=&r" (t)
                        : "r" (location)
                        : "cc", "memory");

   return t;
#else
   return 1 + __fetch_and_add( (volatile int*)location, 1);
#endif  /* !defined(__IBMC__) */
}

#define PASTIX_ATOMIC_HAS_ATOMIC_DEC_32B
static inline int32_t pastix_atomic_dec_32b( volatile int32_t *location )
{
#if !defined(__IBMC__)
   int32_t t;

   __asm__ __volatile__(
                        "1:   lwarx   %0, 0,%1     \n\t"
                        "     addic   %0,%0,-1     \n\t"
                        "     stwcx.  %0,0,%1      \n\t"
                        "     bne-    1b           \n\t"
                        : "=&r" (t)
                        : "r" (location)
                        : "cc", "memory");

   return t;
#else
   return __fetch_and_add( (volatile int*)location, -1) - 1;
#endif  /* !defined(__IBMC__) */
}

#define PASTIX_ATOMIC_HAS_ATOMIC_ADD_32B
static inline int32_t pastix_atomic_add_32b( volatile int32_t *location, int32_t i )
{
#if !defined(__IBMC__)
   int32_t t;

   __asm__ __volatile__(
                        "1:   lwarx   %0, 0, %1    \n\t"
                        "     addic   %0, %0,%2    \n\t"
                        "     stwcx.  %0, 0, %1    \n\t"
                        "     bne-    1b           \n\t"
                        : "=&r" (t)
                        : "r" (location), "r" (i)
                        : "cc", "memory");

   return t;
#else
   return i + __fetch_and_add( (volatile int*)location, i);
#endif  /* !defined(__IBMC__) */
}

#define PASTIX_ATOMIC_HAS_ATOMIC_SUB_32B
static inline int32_t pastix_atomic_sub_32b( volatile int32_t *location, int32_t i )
{
#if !defined(__IBMC__)
   int32_t t;

   __asm__ __volatile__(
                        "1:   lwarx   %0, 0,%1     \n\t"
                        "     subf    %0,%0,%2     \n\t"
                        "     stwcx.  %0,0,%1      \n\t"
                        "     bne-    1b           \n\t"
                        : "=&r" (t)
                        : "r" (location), "r" (i)
                        : "cc", "memory");

   return t;
#else
   return __fetch_and_add( (volatile int*)location, i) - i;
#endif  /* !defined(__IBMC__) */
}

