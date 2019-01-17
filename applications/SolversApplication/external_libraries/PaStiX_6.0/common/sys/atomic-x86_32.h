/*
 * Copyright (c) 2009-2012 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 */

static inline void pastix_mfence(void)
{
    __asm__ __volatile__ ("mfence\n\t":::"memory");
}

static inline int pastix_atomic_cas_32b(volatile uint32_t* location,
                                       uint32_t old_value,
                                       uint32_t new_value)
{
    unsigned char ret;
    __asm__ __volatile__ (
                          "lock; cmpxchgl %3,%4   \n\t"
                          "sete     %0      \n\t"
                          : "=qm" (ret), "=a" (old_value), "=m" (*location)
                          : "q"(new_value), "m"(*location), "1"(old_value)
                          : "memory", "cc");

    return (int)ret;
}

static inline int pastix_atomic_bor_32b(volatile uint32_t* location,
                                       uint32_t value)
{
    uint32_t old_value;

    do {
        old_value = *location;
    } while( !pastix_atomic_cas_32b(location, old_value, (old_value|value) ));
    return old_value | value;
}

static inline int pastix_atomic_band_32b(volatile uint32_t* location,
                                        uint32_t value)
{
    uint32_t old_value;

    do {
        old_value = *location;
    } while( !pastix_atomic_cas_32b(location, old_value, (old_value&value) ));
    return old_value & value;
}

#define ll_low(x)	*(((unsigned int *)&(x)) + 0)
#define ll_high(x)	*(((unsigned int *)&(x)) + 1)

static inline int pastix_atomic_cas_64b(volatile uint64_t* location,
                                       uint64_t old_value,
                                       uint64_t new_value)
{
   /*
    * Compare EDX:EAX with m64. If equal, set ZF and load ECX:EBX into
    * m64. Else, clear ZF and load m64 into EDX:EAX.
    */
    unsigned char ret;

    __asm__ __volatile__(
                    "push %%ebx            \n\t"
                    "movl %3, %%ebx        \n\t"
                    "lock cmpxchg8b (%4)  \n\t"
                    "sete %0               \n\t"
                    "pop %%ebx             \n\t"
                    : "=qm"(ret),"=a"(ll_low(old_value)), "=d"(ll_high(old_value))
                    : "D"(location), "1"(ll_low(old_value)), "2"(ll_high(old_value)),
                      "r"(ll_low(new_value)), "c"(ll_high(new_value))
                    : "cc", "memory", "ebx");
    return (int) ret;
}

#define PASTIX_ATOMIC_HAS_ATOMIC_INC_32B
static inline int32_t pastix_atomic_inc_32b(volatile int32_t *location)
{
    __asm__ __volatile__ (
                          "lock; incl %0\n"
                          : "+m" (*(location)));
    return (*location);
}

#define PASTIX_ATOMIC_HAS_ATOMIC_DEC_32B
static inline int32_t pastix_atomic_dec_32b(volatile int32_t *location)
{
    __asm__ __volatile__ (
                          "lock; decl %0\n"
                          : "+m" (*(location)));
    return (*location);
}

#define PASTIX_ATOMIC_HAS_ATOMIC_ADD_32B
static inline int32_t pastix_atomic_add_32(volatile int32_t* v, int32_t i)
{
    int ret = i;
   __asm__ __volatile__(
                        "lock; xaddl %1,%0"
                        :"=m" (*v), "+r" (ret)
                        :"m" (*v)
                        :"memory", "cc");
   return (ret+i);
}

#define PASTIX_ATOMIC_HAS_ATOMIC_SUB_32B
static inline int32_t pastix_atomic_sub_32(volatile int32_t* v, int32_t i)
{
    int ret = -i;
   __asm__ __volatile__(
                        "lock; xaddl %1,%0"
                        :"=m" (*v), "+r" (ret)
                        :"m" (*v)
                        :"memory", "cc");
   return (ret-i);
}

