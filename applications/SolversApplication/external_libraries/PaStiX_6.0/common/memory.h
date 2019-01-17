/**
 *
 * @file memory.h
 *
 * @copyright 1998-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * PaStiX memory tracking function.
 *
 * @version 6.0.1
 * @author Francois Pellegrini
 * @author Xavier Lacoste
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 */
/*
 *  File: memory.h
 *
 *  Part of a parallel direct block solver.
 *
 *  These lines are the common data
 *  declarations for all modules.
 *
 *  Authors:
 *    Mathieu  Faverge    - faverge@labri.fr
 *    Xavier   LACOSTE    - lacoste@labri.fr
 *    Pierre   RAMET      - ramet@labri.fr
 *
 *  Dates:
 *    Version 0.0 - from 08 may 1998
 *                  to   08 jan 2001
 *    Version 1.0 - from 06 jun 2002
 *                  to   06 jun 2002
 */
#ifndef _memory_h_
#define _memory_h_

/*
 * Function: pastix_protected_malloc
 *
 * PowerPC architectures don't support malloc(0). This function
 * prints a warning when it happens to avoid segfault.
 */
static inline void *pastix_malloc_func( size_t size,
                                        char *filename,
                                        int line )
{
    if (size > 0) {
	return malloc(size);
    }
    else {
	fprintf(stderr, "Pb Alloc 0 %s:%d\n", filename, line);
	return (void *)NULL;
    }
}

#if defined(PASTIX_ARCH_PPC)
#  define memAlloc(size) pastix_malloc_func(size, __FILE__, __LINE__)
#else
#  define memAlloc(size) malloc(size)
#endif

#define memFree(ptr) free((void*)(ptr))
#define memFree_null(ptr) do			\
	{					\
	    memFree( ptr );			\
	    (ptr) = NULL;			\
	} while(0)

#define MALLOC_INTERN(ptr, size, type)                          \
    do {                                                        \
        ptr = (type*)memAlloc((size) * sizeof(type)) ;          \
    } while(0)

#define MALLOC_EXTERN(ptr, size, type)		\
    ptr = (type*)malloc((size) * sizeof(type))

#define MALLOC_ERROR( _str_ )                                           \
    {                                                                   \
        fprintf(stderr, "%s allocation (line=%d,file=%s)\n",(_str_),__LINE__,__FILE__); \
        exit(-1);                                                       \
    }

/*
 * Macro: MALLOC_INTOREXTERN
 *
 * Choose between <MALLOC_INTERN> and <MALLOC_EXTERN>
 * following flag_int.
 *
 * Parameters:
 *   ptr      - address where to allocate.
 *   size     - Number of elements to allocate.
 *   types    - Type of the elements to allocate.
 *   flag_int - 1 for internal allocation, 0 for external.
 */
#define MALLOC_INTOREXTERN(ptr, size, type, flag_int) \
  do {                                                \
    if (flag_int == 1)                          \
      {                                               \
        MALLOC_INTERN(ptr, size, type);               \
      }                                               \
    else                                              \
      {                                               \
        MALLOC_EXTERN(ptr, size, type);               \
      }                                               \
  } while (0)

#define FREE_NULL_INTOREXT(ptr, flag_int)         \
  do {                                            \
    if (flag_int == 1)                      \
      {                                           \
        memFree_null(ptr);                        \
      }                                           \
    else                                          \
      {                                           \
        free(ptr);                                \
        ptr = NULL;                               \
      }                                           \
  } while (0)

#define memRealloc realloc

#endif /* _memory_h_ */
