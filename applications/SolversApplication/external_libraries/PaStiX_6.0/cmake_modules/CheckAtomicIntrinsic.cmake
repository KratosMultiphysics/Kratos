#
# Check if the compiler supports __sync_bool_compare_and_swap.
#
if(NOT HAVE_COMPARE_AND_SWAP_32 AND NOT HAVE_COMPARE_AND_SWAP_64 AND NOT HAVE_COMPARE_AND_SWAP_128)

include(CheckCSourceCompiles)

# Gcc style atomics?
CHECK_C_SOURCE_COMPILES("
      #include <stdint.h>
      int main( int argc, char** argv) {
         int32_t where = 0;
         if (!__sync_bool_compare_and_swap(&where, 0, 1))
            return -1;
         return 0;
      }
      " HAVE_ATOMIC_GCC_32_BUILTINS)
if( HAVE_ATOMIC_GCC_32_BUILTINS )
  CHECK_C_SOURCE_COMPILES("
        #include <stdint.h>
        int main( int argc, char** argv) {
           int64_t where = 0;
           if (!__sync_bool_compare_and_swap(&where, 0, 1))
              return -1;
           return 0;
        }
        " HAVE_ATOMIC_GCC_64_BUILTINS)
endif( HAVE_ATOMIC_GCC_32_BUILTINS )
if( HAVE_ATOMIC_GCC_64_BUILTINS )
  CHECK_C_SOURCE_COMPILES("
        #include <stdint.h>
        int main( int argc, char** argv ) {
            __int128_t where = 0;
            if( !__sync_bool_compare_and_swap(&where, 0, 1))
                return -1;
            return 0;
        }
        " HAVE_ATOMIC_GCC_128_BUILTINS)
  if( NOT HAVE_ATOMIC_GCC_128_BUILTINS ) # try again with -mcx16
    include(CMakePushCheckState)
    CMAKE_PUSH_CHECK_STATE()
    SET( CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -mcx16" )
    UNSET( HAVE_ATOMIC_GCC_128_BUILTINS CACHE )
    CHECK_C_SOURCE_COMPILES("
        #include <stdint.h>
        int main( int argc, char** argv ) {
            __int128_t where = 0;
            if( !__sync_bool_compare_and_swap(&where, 0, 1))
                return -1;
            return 0;
        }
        " HAVE_ATOMIC_GCC_128_BUILTINS)
    CMAKE_POP_CHECK_STATE()
    if( HAVE_ATOMIC_GCC_128_BUILTINS )
      SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcx16" CACHE STRING "" FORCE)
    endif( HAVE_ATOMIC_GCC_128_BUILTINS )
  endif( NOT HAVE_ATOMIC_GCC_128_BUILTINS )
endif( HAVE_ATOMIC_GCC_64_BUILTINS )

# Xlc style atomics?
CHECK_C_SOURCE_COMPILES("
      #include <stdint.h>

      int main( int argc, char** argv)
      {
         int32_t where = 0, old = where;

         if (!__compare_and_swap(&where, &old, 1))
            return -1;

         return 0;
      }
      " HAVE_ATOMIC_XLC_32_BUILTINS)
if( HAVE_ATOMIC_XLC_32_BUILTINS )
  CHECK_C_SOURCE_COMPILES("
        #include <stdint.h>

        int main( int argc, char** argv)
        {
           long where = 0, old = where;

           if (!__compare_and_swaplp(&where, &old, 1))
              return -1;

           return 0;
        }
        " HAVE_ATOMIC_XLC_64_BUILTINS)
endif( HAVE_ATOMIC_XLC_32_BUILTINS )

# MIPS style atomics?
CHECK_C_SOURCE_COMPILES("
      #include <stdint.h>

      int main(int, const char**)
      {
         uint32_t where  = 0;
         if (!__sync_compare_and_swap(&where, 0, 1))
            return -1;

         return 0;
      }
      " HAVE_ATOMIC_MIPOSPRO_32_BUILTINS)
if( HAVE_ATOMIC_MIPOSPRO_32_BUILTINS )
  CHECK_C_SOURCE_COMPILES("
        #include <stdint.h>

        int main(int, const char**)
        {
           uint64_t where  = 0;
           if (!__sync_compare_and_swap(&where, 0, 1))
              return -1;

           return 0;
        }
        " HAVE_ATOMIC_MIPOSPRO_64_BUILTINS)
endif( HAVE_ATOMIC_MIPOSPRO_32_BUILTINS )

# SUN OS style atomics? 
CHECK_C_SOURCE_COMPILES("
      #include <atomic.h>
      #include <stdint.h>

      int main(int, const char**)
      {
         uint_t where = 0;
         if (0 != atomic_cas_uint(&where, 0, 1))
            return -1;

         return 0;
      }
      " HAVE_ATOMIC_SUN_32)
if( HAVE_ATOMIC_SUN_32 )
  CHECK_C_SOURCE_COMPILES("
        #include <atomic.h>
        #include <stdint.h>

        int main(int, const char**)
        {
           uint64_t where = 0;
           if (0 != atomic_cas_uint(&where, 0, 1))
              return -1;

           return 0;
        }
        " HAVE_ATOMIC_SUN_64)
endif( HAVE_ATOMIC_SUN_32 )

# Apple style atomics?
if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  CHECK_FUNCTION_EXISTS(OSAtomicCompareAndSwap32 HAVE_COMPARE_AND_SWAP_32)
  CHECK_FUNCTION_EXISTS(OSAtomicCompareAndSwap64 HAVE_COMPARE_AND_SWAP_64)
endif(CMAKE_SYSTEM_NAME MATCHES "Darwin")

if( HAVE_ATOMIC_SUN_32 OR HAVE_ATOMIC_MIPOSPRO_32_BUILTINS OR HAVE_ATOMIC_GCC_32_BUILTINS )
  set( HAVE_COMPARE_AND_SWAP_32 1 CACHE INTERNAL "Atomic operation on 32 bits are supported")
endif( HAVE_ATOMIC_SUN_32 OR HAVE_ATOMIC_MIPOSPRO_32_BUILTINS OR HAVE_ATOMIC_GCC_32_BUILTINS )

if( HAVE_ATOMIC_SUN_64 OR HAVE_ATOMIC_MIPOSPRO_64_BUILTINS OR HAVE_ATOMIC_GCC_64_BUILTINS )
  set( HAVE_COMPARE_AND_SWAP_64 1 CACHE INTERNAL "Atomic operation on 64 bits are supported")
endif( HAVE_ATOMIC_SUN_64 OR HAVE_ATOMIC_MIPOSPRO_64_BUILTINS OR HAVE_ATOMIC_GCC_64_BUILTINS )

if( HAVE_ATOMIC_GCC_128_BUILTINS )
  set( HAVE_COMPARE_AND_SWAP_128 1 CACHE INTERNAL "Atomic operation on 128 bits are supported")
endif( HAVE_ATOMIC_GCC_128_BUILTINS )

if( HAVE_COMPARE_AND_SWAP_32 )
  message( STATUS "\t support for 32 bits atomics - found")
endif( HAVE_COMPARE_AND_SWAP_32 )

if( HAVE_COMPARE_AND_SWAP_64 )
  message( STATUS "\t support for 64 bits atomics - found")
endif( HAVE_COMPARE_AND_SWAP_64 )

if( HAVE_COMPARE_AND_SWAP_128 )
  message( STATUS "\t support for 128 bits atomics - found")
endif( HAVE_COMPARE_AND_SWAP_128 )

if( CMAKE_SIZEOF_VOID_P MATCHES "8" )
  if( HAVE_COMPARE_AND_SWAP_32 AND NOT HAVE_COMPARE_AND_SWAP_64 )
    message( FATAL_ERROR "64 bits OS with support for 32 bits atomics but without support for 64 bits atomics")
  endif( HAVE_COMPARE_AND_SWAP_32 AND NOT HAVE_COMPARE_AND_SWAP_64 )
  if( NOT HAVE_COMPARE_AND_SWAP_128 )
    message( WARNING "128 bit atomics not found but pointers are 64 bits. Some list operations will not be optimized")
  endif( NOT HAVE_COMPARE_AND_SWAP_128 )
endif( CMAKE_SIZEOF_VOID_P MATCHES "8" )

endif(NOT HAVE_COMPARE_AND_SWAP_32 AND NOT HAVE_COMPARE_AND_SWAP_64 AND NOT HAVE_COMPARE_AND_SWAP_128)
