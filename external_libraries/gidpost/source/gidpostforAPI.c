/* gidpost */

#include <assert.h>
#include <stdarg.h>

#ifdef HAS_INTEL_FORTRAN
#define f2cFortran
#endif
#include "cfortran/cfortran.h"

// Do not complain about implementing deprecated api
#ifdef _WIN32
// #if defined( _MSC_VER )
// disable deprecated declarations
#pragma warning(disable:4996)
// #endif
#else // _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif // _WIN32

#include "gidpost.h"

#ifdef HAS_INTEL_FORTRAN

// TODO
#ifdef __INTEL_LLVM_COMPILER
// Using Intel IFX compiler

// do the same as HAS_GFORTRAN or HAS_NVIDIA_FORTRAN
//e.g with IFX (is also Intel FORTRAN compiles, but old, must be done in this way or it doesn't work at runtime, at least used in IBER GID_OPENPOSTRESULTFILE is not creating the .post.res file)

#include "gidpostfor.h"

#ifdef fcallsc
#undef fcallsc
#endif
#define fcallsc(UN,LN) append_fcallsc(_,_,UN,LN)
#include "gidpostfor.h"


#else // __INTEL_LLVM_COMPILER
// Using Intel IFORT compiler (classic)

#ifdef fcallsc
#undef fcallsc
#endif
// #define fcallsc(UN,LN) LN
// This works for Intel Fortran Classic:
#define fcallsc(UN,LN) UN

#include "gidpostfor.h"


#endif // __INTEL_LLVM_COMPILER


#elif defined(HAS_GFORTRAN) || defined(HAS_NVIDIA_FORTRAN) // HAS_INTEL_FORTRAN

#include "gidpostfor.h"

#ifdef fcallsc
#undef fcallsc
#endif
#define fcallsc(UN,LN) append_fcallsc(_,_,UN,LN)
#include "gidpostfor.h"

#else // defined(HAS_GFORTRAN) || defined(HAS_NVIDIA_FORTRAN) // HAS_INTEL_FORTRAN

// Do the same as GFortran or
// add the proper ifdef

#include "gidpostfor.h"

#ifdef fcallsc
#undef fcallsc
#endif
#define fcallsc(UN,LN) append_fcallsc(_,_,UN,LN)
#include "gidpostfor.h"

#endif

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif // _WIN32
