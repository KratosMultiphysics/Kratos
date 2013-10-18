/* gidpost 1.8 */

#include <assert.h>
#include <stdarg.h>
#ifdef IFORT
#define f2cFortran
#endif
#include "cfortran/cfortran.h"
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
