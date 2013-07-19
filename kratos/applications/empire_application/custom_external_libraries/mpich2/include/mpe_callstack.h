/*
   (C) 2007 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/

#ifndef _MPE_CALLSTACK_H_
#define _MPE_CALLSTACK_H_

#define HAVE_BACKTRACE 1

#define MPE_CALLSTACK_UNLIMITED   9999
#define MPE_CALLSTACK_MAXLINE     1024

#if defined( HAVE_BACKTRACE )

#include <execinfo.h>

#define MPE_CALLSTACK_MAXDEPTH     128

typedef struct MPE_CallStack {
     void    *buffer[ MPE_CALLSTACK_MAXDEPTH ];
     size_t   depth;
     FILE    *pipefile;
     char     line_buf[ MPE_CALLSTACK_MAXLINE ];
} MPE_CallStack_t;

#define MPE_CallStack_init(cstk) \
        do { \
            (cstk)->depth = backtrace( (cstk)->buffer, \
                                       MPE_CALLSTACK_MAXDEPTH ); \
            (cstk)->pipefile = NULL; \
        } while (0)

#define MPE_CallStack_print(cstk,fd) \
        do { \
            backtrace_symbols_fd( (cstk)->buffer, (cstk)->depth, (fd) ); \
        } while (0)

#elif defined( HAVE_PRINTSTACK )

#include <ucontext.h>

typedef struct MPE_CallStack {
     FILE    *pipefile;
     char     line_buf[ MPE_CALLSTACK_MAXLINE ];
} MPE_CallStack_t;

#define MPE_CallStack_init(cstk)

#define MPE_CallStack_print(cstk,fd) \
        do { \
            printstack( (fd) ); \
        } while (0)

#else




#if defined( HAVE_STRING_H )
#include <string.h>
#endif
#if defined( HAVE_UNISTD_H )
#include <unistd.h>
#endif

typedef struct MPE_CallStack {
     FILE    *pipefile;
     char     line_buf[ MPE_CALLSTACK_MAXLINE ];
} MPE_CallStack_t;

#define MPE_CallStack_init(cstk)

#define MPE_CallStack_print(cstk,fd) \
        do { \
            char msg[] = "No callstack trace available.\n"; \
            write( fd, msg, strlen(msg) ); \
        } while (0)

#endif

/* void MPE_CallStack_init( MPE_CallStack_t *cstk ) */
/* void MPE_CallStack_print( MPE_CallStack_t *cstk, int fd ) */
void        MPE_CallStack_iteratorInit( MPE_CallStack_t *cstk );
int         MPE_CallStack_iteratorHasMore( MPE_CallStack_t *cstk );
const char* MPE_CallStack_iteratorFetchNext( MPE_CallStack_t *cstk );

void        MPE_CallStack_fancyprint( MPE_CallStack_t *cstk, int fd,
                                      const char      *line_prefix,
                                            int        to_print_stack_index,
                                            int        max_printed_frames );

#endif
