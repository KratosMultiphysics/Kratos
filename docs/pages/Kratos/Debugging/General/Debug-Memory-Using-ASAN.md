---
title: Debug Memory Using ASAN
keywords: 
tags: [Debug-Memory-Using-ASAN.md]
sidebar: kratos_debugging
summary: 
---

_**Note**: Please note this page assumes advance knowledge of compilers._

## Debug Memory Using ASAN

It is possible to compile Kratos with different sanitizers which can point to different memory problems that otherwise are hard to locate. As this conflicts with some basic compilations options enabled by default in Kratos cannot be directly added as a separate build mode. 
This process has been tested with clang, but gcc should be suitable for it as well.
The steps to compile are detailed here:

### CMake Changes:

- Remove `-Wl,--no-undefined` from all compilation targets. This must be done in the root `CMakeLists.txt`
- Add `link_libraries(-fsanitize=address -shared-libasan)` before declaring any target_link_library. Ideally you want to do that before line 200 of the root `CMakeLists.txt`.

Example of a possible diff:

```CMake
# Set compiler flags
if(${CMAKE_COMPILER_IS_GNUCXX})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops -Wall -std=c++11")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override")
  endif()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funroll-loops -Wall")
  # set(CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")
  message(STATUS "additional default options were set for gcc")
  message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
  message(STATUS "CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
endif(${CMAKE_COMPILER_IS_GNUCXX})

message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  	  if(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -funroll-loops -Wall -Wno-unused-local-typedef -Wno-unknown-pragmas  ")
		set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -funroll-loops -Wall -Wno-unknown-pragmas  ")
		# set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
		message("additional default options were set for clang compiler")
		message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
		message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
	  endif(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
  else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	  if(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -funroll-loops -Wall -Wno-unused-local-typedef -Wno-unknown-pragmas  ")
		set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops -Wall -Wno-unknown-pragmas  ")
		# set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
		message("additional default options were set for clang compiler")
		message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
		message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
	  endif(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
	  if(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
		if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0)
		  # pybind requires min. version 17, but we need at least 18:
		  message( "DEPRECATED: detected compiler as Intel " ${CMAKE_CXX_COMPILER_VERSION} )
		  message( FATAL_ERROR "Please use Version 18 or greater")
		endif()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC  -funroll-loops -lpthread -wd654 -wd10010  ")
		set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops -lpthread -wd654 -wd10010  ")
		# set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
		message("additional default options were set for intel compiler")
		message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
		message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
	  endif(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
else(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  if(${CMAKE_C_COMPILER} MATCHES "icc.*$")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC  -funroll-loops  ") #-lpthread -wd654 -wd10010  ")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops  ") #-lpthread -wd654 -wd10010  ")
    # set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
    message("additional default options were set for intel compiler")
    message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
    message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
  endif(${CMAKE_C_COMPILER} MATCHES "icc.*$")
endif(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

link_libraries(-fsanitize=address -shared-libasan)
```

### Configure.sh Changes

You need to add the following compiler flags:
- -fsanitize=address: Enables the memory sanitizer
- -fno-omit-frame-pointer: Allows you to obtain full stack-traces
- -fsanitize-recover=address: Allows the sanitizer to recover from errors (allows execution to continue after the first detection)

### Symbolizer
In order for the address sanitizer to give you code names instead of addresses you will need a symbolyzer. You can find `llvm-symbolyzer` as part of the `llvm` package or as stand-alone program in some package managers. Once installed:

```
export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer
```

### Launch the code
Before running you must set some options:

#### ASAN_OPTIONS

- suppressions=pybind.sup: Allows yo define a suppression file (similar to valgrind ones). Currently appears to be bugged while launching from python.

- halt_on_error=0/1: If compiled with `-fsanitize-recover=address` allows the code to continue after detecting the first problem.

- fast_unwind_on_malloc=0/1: If set to 1 increases speed at the cost of preventing part of the trace from being collected. It is recommended to set it to 0.

Example:
```
export ASAN_OPTIONS=suppressions=pybind.supp:halt_on_error=0:fast_unwind_on_malloc=0
```

It is also necessary to preload the asan library if you are running from python. In order to do so:

For serial
```
LD_PRELOAD=/usr/lib/clang/11.0.0/lib/linux/libclang_rt.asan-x86_64.so python script.py
```

For MPI:
```
mpirun -x LD_PRELOAD=/usr/lib/clang/11.0.0/lib/linux/libclang_rt.asan-x86_64.so -np 2 --output-filename asan_test python script.py

```


