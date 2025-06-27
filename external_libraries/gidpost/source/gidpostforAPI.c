/* gidpost */

#include <assert.h>
#include <stdarg.h>
#ifdef IFORT
#define f2cFortran
#endif
#include "cfortran/cfortran.h"

// Do not complain about implementing deprecated api
#ifdef WIN32
// #if defined( _MSC_VER )
// disable deprecated declarations
#pragma warning(disable:4996)
// #endif
#else // WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif // WIN32

#include "gidpost.h"

#ifdef IFORT

#ifdef fcallsc
#undef fcallsc
#endif
#define fcallsc(UN,LN) LN

#include "gidpostfor.h"

#else

#include "gidpostfor.h"

#ifdef fcallsc
#undef fcallsc
#endif
#define fcallsc(UN,LN) append_fcallsc(_,_,UN,LN)

#include "gidpostfor.h"

#endif

#ifndef WIN32
#pragma GCC diagnostic pop
#endif // WIN32
