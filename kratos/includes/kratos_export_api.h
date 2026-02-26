//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig Pina
//

// System includes

// External includes

// Project includes

#pragma once

#undef KRATOS_API_EXPORT
#undef KRATOS_API_IMPORT
#if defined(_WIN32) || defined(_WIN64)
    #if defined(__MINGW32__) || defined(__MINGW64__) || defined(__INTEL_LLVM_COMPILER )
        #define KRATOS_API_EXPORT __attribute__((visibility("default")))
        #define KRATOS_API_IMPORT __attribute__((visibility("default")))
    #else
        #define KRATOS_API_EXPORT __declspec(dllexport)
        #define KRATOS_API_IMPORT __declspec(dllimport)
    #endif
#else
    #define KRATOS_API_EXPORT __attribute__((visibility("default")))
    #define KRATOS_API_IMPORT __attribute__((visibility("default")))
#endif

// Fixes MSVC not expanding __VA_ARGS__ as defined in the C99 standard
#define KRATOS_EXPAND(A) A

// Expands the API call to either import or export based on the
// number of the arguments in API() call.
#define KRATOS_API_CALL(x,T1,T2,T3,...) T3

// If KRATOS_API_NO_DLL is defined ignore the DLL API
#ifndef KRATOS_API_NO_DLL
    // Intel dpcpp and msvc treatment of __VA_ARGS__ is non-conforming
    #if (defined(_WIN32) || defined(_WIN64)) && defined(__INTEL_LLVM_COMPILER )
        #define KRATOS_API(...) \
            KRATOS_EXPAND(KRATOS_API_CALL(,__VA_ARGS__,KRATOS_API_EXPORT,KRATOS_API_IMPORT))
    #else
        #define KRATOS_API(...) \
            KRATOS_EXPAND(KRATOS_API_CALL(,##__VA_ARGS__,KRATOS_API_EXPORT,KRATOS_API_IMPORT))
    #endif
    #define KRATOS_NO_EXPORT(...)
#else
    #define KRATOS_API(...)
    #define KRATOS_NO_EXPORT(...)
#endif

// Conditionally declare explicit template instances, since explicit instantiation does not play nice with dllexport
#undef KRATOS_API_EXTERN
#if defined(_WIN32) || defined(_WIN64)
    #if defined(__MINGW32__) || defined(__MINGW64__) || defined(__INTEL_LLVM_COMPILER )
        #define KRATOS_API_EXTERN extern
    #else
        #define KRATOS_API_EXTERN
    #endif
#else
    #define KRATOS_API_EXTERN extern
#endif
