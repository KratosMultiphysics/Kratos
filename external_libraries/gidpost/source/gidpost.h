/*
 *  gidpost.h --
 *
 *    This file declare a C interface for generating postprocessing
 *    results in the 'New postprocess format' of GiD. The interface is
 *    composed by a set of functions according to the specification of
 *    the GiD postprocess file format from version 7.2. Althought the
 *    specification only describe the ascii format, with this
 *    interface is possible to write in ascii, ascii zipped and binary
 *    format. Binary files are always zipped. All functions returns an
 *    integer indicating fail with a non zero (!0) value or success
 *    with a zero value (0).
 *
 */

#ifndef __GIDPOST__
#define __GIDPOST__

/* if this is undefined, it is not necessary to link with library hdf5 and format GiD_PostHDF5 gets deactivated
   this macro should be defined in the project configuration
*/
/*
#define ENABLE_HDF5
*/

#ifdef HAVE_GIDPOST_CONFIG_H
#include "gidpost_config.h"
#endif

#define GP_VERSION_MAJOR 2
#define GP_VERSION_MINOR 13

/* build GIDPOST_VERSION by stringify-ing GP_VERSION_MAJOR and GP_VERSION_MINOR */
#define __GP_V_STR( str)     __GP_V_STR_( str)
#define __GP_V_STR_( str)     #str
#define GIDPOST_VERSION     "" __GP_V_STR( GP_VERSION_MAJOR) "." __GP_V_STR( GP_VERSION_MINOR) ""

#define GP_CONST const

#define GP_UNKNOWN    ( ( float)-3.40282346638528860e+38)

/*
#if defined(_MSC_VER) && defined(GIDPOST_SHARED)
# ifndef _WIN32
#  define _WIN32
# endif
# if defined(GIDPOST_EXPORTS)
#  define GIDPOST_API __declspec(dllexport)
# else
#  define GIDPOST_API __declspec(dllimport)
# endif
#else
# define GIDPOST_API
#endif
*/

#if defined (_WIN32)
// in Windows (VisualStudio), deprecated declarations are of the form:
// __declspec(deprecated) void func1() {}
  #if defined (GIDPOST_SHARED)
    #if defined (GIDPOST_EXPORTS)
      #define  GIDPOST_API               __declspec( dllexport)
      #define  GIDPOST_API_DEPRECATED    __declspec( deprecated, dllexport)
    #else
      #define  GIDPOST_API               __declspec( dllimport)
      #define  GIDPOST_API_DEPRECATED    __declspec( deprecated, dllimport)
    #endif /* defined (GIDPOST_EXPORTS) */
  #else // GIDPOST_SHARED
    #define  GIDPOST_API
    #define  GIDPOST_API_DEPRECATED    __declspec( deprecated)
  #endif // GIDPOST_SHARED
  #define GCC_GIDPOST_API_DEPRECATED
#else /* defined (_WIN32)*/
  // in linux / macOS, deprecated declarations are of the form:
  // void __attribute__(( deprecated)) func1() {}
  #define GIDPOST_API
  #define GIDPOST_API_DEPRECATED
  #define GCC_GIDPOST_API_DEPRECATED   __attribute__(( deprecated))
#endif

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif
  
/* ################################################################################
 * #    command to review contents of the file in GiD_PostHDF5 format:
 * #           h5dump --string FILE.flavia.res > FILE.txt
 * #
 * ################################################################################
*/

#include "gidpost_types.h"
#include "gidpost_functions.h"
#include "gidpost_functions_deprecated.h"

#ifdef __cplusplus
}
#endif

#endif
