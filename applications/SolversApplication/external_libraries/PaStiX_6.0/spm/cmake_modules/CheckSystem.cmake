# cmake modules setup
cmake_minimum_required (VERSION 2.8.7)
include (CMakeDetermineSystem)
include (CheckCCompilerFlag)
include (CheckFunctionExists)
include (CheckSymbolExists)
include (CheckIncludeFiles)
include (CMakePushCheckState)

#
# Fix the building system for 32 or 64 bits.
#
# On MAC OS X there is an easy solution, by setting the
# CMAKE_OSX_ARCHITECTURES to a subset of the following values:
# ppc;ppc64;i386;x86_64.
# On Linux this is a little bit tricky. We have to check that the
# compiler supports the -m32/-m64 flags as well as the linker.
# Once this issue is resolved the CMAKE_C_FLAGS and CMAKE_C_LDFLAGS
# have to be updated accordingly.
#
# TODO: Same trick for the Fortran compiler...
#       no idea how to correctly detect if the required/optional
#          libraries are in the correct format.
#
if (BUILD_64bits)
  if( _match_xlc)
    set( ARCH_BUILD "-q64" )
  else (_match_xlc)
    if( ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "sparc64fx" )
      set ( ARCH_BUILD " " )
    else()
      set( ARCH_BUILD "-m64" )
    endif()
  endif(_match_xlc)
else (BUILD_64bits)
  if( _match_xlc)
    set( ARCH_BUILD "-q32" )
  else (_match_xlc)
    set( ARCH_BUILD "-m32" )
  endif(_match_xlc)
endif (BUILD_64bits)

check_c_compiler_flag( ${ARCH_BUILD} C_M32or64 )
if( C_M32or64 )
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ARCH_BUILD}" )
  set( CMAKE_C_LDFLAGS "${CMAKE_C_LDFLAGS} ${ARCH_BUILD}" )
  set( LOCAL_FORTRAN_LINK_FLAGS "${LOCAL_FORTRAN_LINK_FLAGS} ${ARCH_BUILD}" )
endif( C_M32or64 )

# Set warnings for debug builds
check_c_compiler_flag( "-Wall" HAVE_WALL )
if( HAVE_WALL )
    set( C_WFLAGS "${C_WFLAGS} -Wall" )
endif( HAVE_WALL )
check_c_compiler_flag( "-Wextra" HAVE_WEXTRA )
if( HAVE_WEXTRA )
    set( C_WFLAGS "${C_WFLAGS} -Wextra" )
endif( HAVE_WEXTRA )

#
# flags for Intel icc
#
string(REGEX MATCH ".*icc$" _match_icc ${CMAKE_C_COMPILER})
if(_match_icc)
  # Silence annoying warnings
  check_c_compiler_flag( "-wd424" HAVE_WD )
  if( HAVE_WD )
    # 424: checks for duplicate ";"
    # 981: every volatile triggers a "unspecified evaluation order", obnoxious
    #      but might be useful for some debugging sessions.
    # 1419: warning about extern functions being declared in .c
    #       files
    # 1572: cuda compares floats with 0.0f.
    set( C_WFLAGS "${C_WFLAGS} -wd424 -wd981 -wd1419 -wd1572" )
  endif( HAVE_WD )
endif(_match_icc)

set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_WFLAGS}" )

# add gdb symbols in debug and relwithdebinfo
check_c_compiler_flag( "-g3" HAVE_G3 )
if( HAVE_G3 )
    set( CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g3" )
    set( CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -g3" )
else()
    set( CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0" )
endif( HAVE_G3 )

#
# Remove all duplicates from the CFLAGS.
#
set(TMP_LIST ${CMAKE_C_FLAGS})
separate_arguments(TMP_LIST)
list(REMOVE_DUPLICATES TMP_LIST)
set(CMAKE_C_FLAGS "")
foreach( ITEM ${TMP_LIST})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ITEM}")
endforeach()

# Check for attribute fallthrough
# -------------------------------
check_c_source_compiles("
#include <stdarg.h>
int main(void) {
  int a = 2;
  switch( a ){
  case 0:
     __attribute__((fallthrough));
  default:
      ;
  }
  return 0;
}
"
  HAVE_FALLTHROUGH
  )

# Check for Thread library
# ------------------------
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads)
if(Threads_FOUND)
  CMAKE_PUSH_CHECK_STATE()
  list( APPEND CMAKE_REQUIRED_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}" )
  check_function_exists(pthread_create HAVE_PTHREAD)
  CMAKE_POP_CHECK_STATE()
  if(HAVE_PTHREAD)
    list(APPEND EXTRA_LIBS ${CMAKE_THREAD_LIBS_INIT})
  endif(HAVE_PTHREAD)
endif(Threads_FOUND)

check_function_exists(sched_setaffinity HAVE_SCHED_SETAFFINITY)
if( NOT HAVE_SCHED_SETAFFINITY )
  check_library_exists(rt sched_setaffinity "" HAVE_SCHED_SETAFFINITY)
endif( NOT HAVE_SCHED_SETAFFINITY )

# timeval, timespec, realtime clocks, etc
include(CheckStructHasMember)
check_struct_has_member("struct timespec" tv_nsec time.h HAVE_TIMESPEC_TV_NSEC)
if( NOT HAVE_TIMESPEC_TV_NSEC )
  add_definitions(-D_GNU_SOURCE)
  check_struct_has_member("struct timespec" tv_nsec time.h HAVE_TIMESPEC_TV_NSEC)
endif( NOT HAVE_TIMESPEC_TV_NSEC )
check_library_exists(rt clock_gettime "" HAVE_CLOCK_GETTIME)
if( HAVE_CLOCK_GETTIME )
  list(APPEND EXTRA_LIBS rt)
endif( HAVE_CLOCK_GETTIME )

# stdlib, stdio, string, getopt, etc
check_include_files(stdarg.h HAVE_STDARG_H)
# va_copy is special as it is not required to be a function.
if (HAVE_STDARG_H)
  check_c_source_compiles("
      #include <stdarg.h>
      int main(void) {
          va_list a, b;
          va_copy(a, b);
          return 0;
      }"
      HAVE_VA_COPY
      )

  if (NOT HAVE_VA_COPY)
    check_c_source_compiles("
        #include <stdarg.h>
        int main(void) {
            va_list a, b;
            __va_copy(a, b);
            return 0;
        }"
        HAVE_UNDERSCORE_VA_COPY
        )
  endif (NOT HAVE_VA_COPY)
endif (HAVE_STDARG_H)

check_function_exists(asprintf HAVE_ASPRINTF)
check_function_exists(vasprintf HAVE_VASPRINTF)
check_include_files(getopt.h HAVE_GETOPT_H)
check_include_files(unistd.h HAVE_UNISTD_H)
check_function_exists(getopt_long HAVE_GETOPT_LONG)
check_include_files(errno.h HAVE_ERRNO_H)
check_include_files(stddef.h HAVE_STDDEF_H)
check_include_files(stdbool.h HAVE_STDBOOL_H)
check_function_exists(getrusage HAVE_GETRUSAGE)
check_symbol_exists(RUSAGE_THREAD sys/resource.h HAVE_RUSAGE_THREAD)
check_include_files(limits.h HAVE_LIMITS_H)
check_include_files(string.h HAVE_STRING_H)
check_include_files(libgen.h HAVE_GEN_H)
check_include_files(complex.h HAVE_COMPLEX_H)
check_include_files(sys/param.h HAVE_SYS_PARAM_H)
check_include_files(sys/types.h HAVE_SYS_TYPES_H)
check_include_files(syslog.h HAVE_SYSLOG_H)

#
# Fortran tricks
#
IF (CMAKE_Fortran_COMPILER_WORKS)
  STRING(REGEX MATCH "Intel" _match_intel ${CMAKE_Fortran_COMPILER_ID})
  IF (_match_intel)
    MESSAGE(STATUS "Add -nofor_main to the Fortran linker.")
    SET(LOCAL_FORTRAN_LINK_FLAGS "${LOCAL_FORTRAN_LINK_FLAGS} -nofor_main")
  ENDIF (_match_intel)

  STRING(REGEX MATCH "PGI$" _match_pgi ${CMAKE_Fortran_COMPILER_ID})
  IF (_match_pgi)
    MESSAGE(STATUS "Add -Mnomain to the Fortran linker.")
    SET(LOCAL_FORTRAN_LINK_FLAGS "${LOCAL_FORTRAN_LINK_FLAGS} -Mnomain -Bstatic")
  ENDIF (_match_pgi)

  STRING(REGEX MATCH ".*xlc$" _match_xlc ${CMAKE_C_COMPILER})
  IF (_match_xlc)
    MESSAGE(ERROR "Please use the thread-safe version of the xlc compiler (xlc_r)")
  ENDIF (_match_xlc)
  STRING(REGEX MATCH "XL" _match_xlc ${CMAKE_C_COMPILER_ID})
  IF (_match_xlc AND BUILD_64bits)
    MESSAGE(STATUS "Add -q64 to the C compiler/linker.")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -q64")
  ENDIF (_match_xlc AND BUILD_64bits)

  STRING(REGEX MATCH ".*xlf$" _match_xlf ${CMAKE_Fortran_COMPILER})
  IF (_match_xlf)
    MESSAGE(ERROR "Please use the thread-safe version of the xlf compiler (xlf_r)")
  ENDIF (_match_xlf)
  STRING(REGEX MATCH "XL$" _match_xlf ${CMAKE_Fortran_COMPILER_ID})
  IF (_match_xlf)
    SET(arch_flags "-q32")
    IF(BUILD_64bits)
      SET(arch_flags "-q64")
    ENDIF(BUILD_64bits)
    MESSAGE(STATUS "Add ${arch_flags} and -nofor_main to the Fortran linker.")
    SET(LOCAL_FORTRAN_LINK_FLAGS "${LOCAL_FORTRAN_LINK_FLAGS} ${arch_flags} -nofor_main")
  ENDIF (_match_xlf)

#
# Even more Fortran tricks.
#
# FFLAGS depend on the compiler

  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    # gfortran
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
    list(APPEND EXTRA_LIBS gfortran)
  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    # ifort
    set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
    string (REPLACE "-i_dynamic" "" CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS}")
  else (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
ENDIF (CMAKE_Fortran_COMPILER_WORKS)
