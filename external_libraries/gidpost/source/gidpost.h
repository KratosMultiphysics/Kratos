/* gidpost 2.11 */
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
#define GP_VERSION_MINOR 11

/* build GIDPOST_VERSION by stringify-ing GP_VERSION_MAJOR and GP_VERSION_MINOR */
#define __GP_V_STR( str)     __GP_V_STR_( str)
#define __GP_V_STR_( str)     #str
#define GIDPOST_VERSION     "" __GP_V_STR( GP_VERSION_MAJOR) "." __GP_V_STR( GP_VERSION_MINOR) ""

#define GP_CONST const

#define GP_UNKNOWN    ( ( float)-3.40282346638528860e+38)

/*
#if defined(_MSC_VER) && defined(GIDPOST_SHARED)
# ifndef WIN32
#  define WIN32
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

// First we define the CPP_DEPRECATED macro for marking deprecated functions.
#if __cplusplus >= 201402L // C++14 ou plus
    #define CPP_DEPRECATED(msg) [[deprecated(msg)]]
#elif defined(_MSC_VER)
    #define CPP_DEPRECATED(msg) __declspec(deprecated(msg))
#elif defined(__GNUC__) || defined(__clang__)
    #define CPP_DEPRECATED(msg) __attribute__((deprecated(msg)))
#else
    #define CPP_DEPRECATED(msg)
#endif

// Now define the GIDPOST_API and GIDPOST_API_DEPRECATED macros
#if defined (WIN32)
  #if defined (GIDPOST_SHARED)
    #if defined (GIDPOST_EXPORTS)
      #define GIDPOST_API __declspec(dllexport)
      #define GIDPOST_API_DEPRECATED CPP_DEPRECATED("This function is deprecated") __declspec(dllexport)
    #else
      #define GIDPOST_API __declspec(dllimport)
      #define GIDPOST_API_DEPRECATED CPP_DEPRECATED("This function is deprecated") __declspec(dllimport)
    #endif
  #else
    #define GIDPOST_API
    #define GIDPOST_API_DEPRECATED CPP_DEPRECATED("This function is deprecated")
  #endif
#else // Non-WIN32
  #if __GNUC__ >= 4 // GCC 4+ supports visibility
    #define GIDPOST_API __attribute__((visibility("default")))
  #else
    #define GIDPOST_API
  #endif
  #define GIDPOST_API_DEPRECATED CPP_DEPRECATED("This function is deprecated") GIDPOST_API // Combine deprecation and visibility if necessary
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
