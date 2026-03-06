//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#pragma once

/* System includes */
#include <stdexcept>
#include <sstream>

/* External includes */

/* Project includes */
#include "includes/kratos_export_api.h"
#include "includes/smart_pointers.h"
#include "includes/exception.h"

// Defining the OS
#if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    #define KRATOS_COMPILED_IN_LINUX
#elif defined(__APPLE__) && defined(__MACH__)
    #define KRATOS_COMPILED_IN_OS
#elif defined(_WIN32) || defined(_WIN64)
    #define KRATOS_COMPILED_IN_WINDOWS
#endif

// Defining the architecture (see https://sourceforge.net/p/predef/wiki/Architectures/)
// Check Windows
#if defined(_WIN32) || defined(_WIN64)
   #if defined(_WIN64)
     #define KRATOS_ENV64BIT
   #else
     #define KRATOS_ENV32BIT
     #error 32 bit system are not supported anymore. Please consider a 64 bits system
  #endif
#else // It is POSIX (Linux, MacOSX, BSD...)
  #if defined(__x86_64__) || defined(__ppc64__) || defined(__aarch64__)
    #define KRATOS_ENV64BIT
  #else // This includes __arm__ and __x86__
    #define KRATOS_ENV32BIT
     #error 32 bit system are not supported anymore. Please consider a 64 bits system
  #endif
#endif

//-----------------------------------------------------------------
//
// Warnings
//
//-----------------------------------------------------------------

#if defined(_MSC_VER)
#  pragma warning(disable: 4244 4267)
#endif

//-----------------------------------------------------------------
//
// Exceptions
//
//-----------------------------------------------------------------

#define KRATOS_CATCH_AND_THROW(ExceptionType, MoreInfo, Block) \
catch(ExceptionType& e)                                        \
{                                                              \
Block                                                          \
KRATOS_ERROR << e.what();                             \
}

#define KRATOS_THROW_ERROR(ExceptionType, ErrorMessage, MoreInfo)    \
{                                                              \
KRATOS_ERROR << ErrorMessage << MoreInfo << std::endl;          \
}

#define KRATOS_CATCH_WITH_BLOCK(MoreInfo,Block) \
} \
KRATOS_CATCH_AND_THROW(std::overflow_error,MoreInfo,Block)   \
KRATOS_CATCH_AND_THROW(std::underflow_error,MoreInfo,Block)  \
KRATOS_CATCH_AND_THROW(std::range_error,MoreInfo,Block)      \
KRATOS_CATCH_AND_THROW(std::out_of_range,MoreInfo,Block)     \
KRATOS_CATCH_AND_THROW(std::length_error,MoreInfo,Block)     \
KRATOS_CATCH_AND_THROW(std::invalid_argument,MoreInfo,Block) \
KRATOS_CATCH_AND_THROW(std::domain_error,MoreInfo,Block)     \
KRATOS_CATCH_AND_THROW(std::logic_error,MoreInfo,Block)      \
KRATOS_CATCH_AND_THROW(std::runtime_error,MoreInfo,Block)    \
catch(Exception& e) { Block throw Exception(e) << KRATOS_CODE_LOCATION << MoreInfo << std::endl; } \
catch(std::exception& e) { Block KRATOS_THROW_ERROR(std::runtime_error, e.what(), MoreInfo) } \
catch(...) { Block KRATOS_THROW_ERROR(std::runtime_error, "Unknown error", MoreInfo) }

#define KRATOS_CATCH_BLOCK_BEGIN class ExceptionBlock{public: void operator()(void){
#define KRATOS_CATCH_BLOCK_END }} exception_block; exception_block();

#ifndef KRATOS_NO_TRY_CATCH
    #define KRATOS_TRY_IMPL try {
    #define KRATOS_CATCH_IMPL(MoreInfo) KRATOS_CATCH_WITH_BLOCK(MoreInfo,{})
#else
    #define KRATOS_TRY_IMPL {};
    #define KRATOS_CATCH_IMPL(MoreInfo) {};
#endif

#ifndef __SUNPRO_CC
    #define KRATOS_TRY KRATOS_TRY_IMPL
    #define KRATOS_CATCH(MoreInfo) KRATOS_CATCH_IMPL(MoreInfo)
#else
    #define KRATOS_TRY {};
    #define KRATOS_CATCH(MoreInfo) {};
#endif

//-----------------------------------------------------------------
//
// macro argument manipulation
//
//-----------------------------------------------------------------

#undef CAT_
#define CAT_(a, b) a##_##b

#undef STR_
#define STR_(x) #x

#undef STR
#define STR(x) STR_(x)

#undef CAT_STR
#define CAT_STR(a, b) STR(CAT_(a, b))

//-----------------------------------------------------------------
//
// variables
//
//-----------------------------------------------------------------

#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

#undef KRATOS_DEFINE_VARIABLE
#define KRATOS_DEFINE_VARIABLE(type, name) \
    inline constexpr Variable<type> name = Variable<type>(#name);

#undef KRATOS_DEFINE_APPLICATION_VARIABLE
#define KRATOS_DEFINE_APPLICATION_VARIABLE(application, type, name) \
    inline constexpr Variable<type> name = Variable<type>(#name);

#undef KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS
#define KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(name) \
    inline constexpr Variable<Kratos::array_1d<double, 3>> name = Variable<Kratos::array_1d<double, 3>>(#name); \
    inline constexpr Variable<double> CAT_(name, X) = Variable<double>(CAT_STR(name, X), &name, 0); \
    inline constexpr Variable<double> CAT_(name, Y) = Variable<double>(CAT_STR(name, Y), &name, 1); \
    inline constexpr Variable<double> CAT_(name, Z) = Variable<double>(CAT_STR(name, Z), &name, 2);

#undef KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS
#define KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(application, name) \
    inline constexpr Variable<Kratos::array_1d<double, 3>> name = Variable<Kratos::array_1d<double, 3>>(#name); \
    inline constexpr Variable<double> CAT_(name, X) = Variable<double>(CAT_STR(name, X), &name, 0); \
    inline constexpr Variable<double> CAT_(name, Y) = Variable<double>(CAT_STR(name, Y), &name, 1); \
    inline constexpr Variable<double> CAT_(name, Z) = Variable<double>(CAT_STR(name, Z), &name, 2);

#undef KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#define KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
    inline constexpr Variable<Kratos::array_1d<double, 3>> name = Variable<Kratos::array_1d<double, 3>>(#name); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, YY), &name, 1); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 2);

#undef KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS
#define KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS(application, name) \
    inline constexpr Variable<Kratos::array_1d<double, 3>> name = Variable<Kratos::array_1d<double, 3>>(#name); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, YY), &name, 1); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 2);

#undef KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#define KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
    inline constexpr Variable<Kratos::array_1d<double, 6>> name = Variable<Kratos::array_1d<double, 6>>(#name); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, YY), &name, 1); \
    inline constexpr Variable<double> CAT_(name, ZZ) = Variable<double>(CAT_STR(name, ZZ), &name, 2); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 3); \
    inline constexpr Variable<double> CAT_(name, YZ) = Variable<double>(CAT_STR(name, YZ), &name, 4); \
    inline constexpr Variable<double> CAT_(name, XZ) = Variable<double>(CAT_STR(name, XZ), &name, 5);

#undef KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS
#define KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS(application, name) \
    inline constexpr Variable<Kratos::array_1d<double, 6>> name = Variable<Kratos::array_1d<double, 6>>(#name); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, YY), &name, 1); \
    inline constexpr Variable<double> CAT_(name, ZZ) = Variable<double>(CAT_STR(name, ZZ), &name, 2); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 3); \
    inline constexpr Variable<double> CAT_(name, YZ) = Variable<double>(CAT_STR(name, YZ), &name, 4); \
    inline constexpr Variable<double> CAT_(name, XZ) = Variable<double>(CAT_STR(name, XZ), &name, 5);

#undef KRATOS_DEFINE_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#define KRATOS_DEFINE_2D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
    inline constexpr Variable<Kratos::array_1d<double, 4>> name = Variable<Kratos::array_1d<double, 4>>(#name); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, YY), &name, 1); \
    inline constexpr Variable<double> CAT_(name, YX) = Variable<double>(CAT_STR(name, ZZ), &name, 2); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, XY), &name, 3);

#undef KRATOS_DEFINE_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS
#define KRATOS_DEFINE_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS(application, name) \
    inline constexpr Variable<Kratos::array_1d<double, 4>> name = Variable<Kratos::array_1d<double, 4>>(#name); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, YY), &name, 1); \
    inline constexpr Variable<double> CAT_(name, YX) = Variable<double>(CAT_STR(name, ZZ), &name, 2); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, XY), &name, 3);

#undef KRATOS_DEFINE_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#define KRATOS_DEFINE_3D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
    inline constexpr Variable<Kratos::array_1d<double, 9>> name = Variable<Kratos::array_1d<double, 9>>(#name); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 1); \
    inline constexpr Variable<double> CAT_(name, XZ) = Variable<double>(CAT_STR(name, XZ), &name, 2); \
    inline constexpr Variable<double> CAT_(name, YX) = Variable<double>(CAT_STR(name, YX), &name, 3); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, YY), &name, 4); \
    inline constexpr Variable<double> CAT_(name, YZ) = Variable<double>(CAT_STR(name, YZ), &name, 5); \
    inline constexpr Variable<double> CAT_(name, ZX) = Variable<double>(CAT_STR(name, ZX), &name, 6); \
    inline constexpr Variable<double> CAT_(name, ZY) = Variable<double>(CAT_STR(name, ZY), &name, 7); \
    inline constexpr Variable<double> CAT_(name, ZZ) = Variable<double>(CAT_STR(name, ZZ), &name, 8);

#undef KRATOS_DEFINE_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS
#define KRATOS_DEFINE_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS(application, name) \
    inline constexpr Variable<Kratos::array_1d<double, 9>> name = Variable<Kratos::array_1d<double, 9>>(#name); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 1); \
    inline constexpr Variable<double> CAT_(name, XZ) = Variable<double>(CAT_STR(name, XZ), &name, 2); \
    inline constexpr Variable<double> CAT_(name, YX) = Variable<double>(CAT_STR(name, YX), &name, 3); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, YY), &name, 4); \
    inline constexpr Variable<double> CAT_(name, YZ) = Variable<double>(CAT_STR(name, YZ), &name, 5); \
    inline constexpr Variable<double> CAT_(name, ZX) = Variable<double>(CAT_STR(name, ZX), &name, 6); \
    inline constexpr Variable<double> CAT_(name, ZY) = Variable<double>(CAT_STR(name, ZY), &name, 7); \
    inline constexpr Variable<double> CAT_(name, ZZ) = Variable<double>(CAT_STR(name, ZZ), &name, 8);

#undef KRATOS_CREATE_VARIABLE
#define KRATOS_CREATE_VARIABLE(type, name) ;

#undef KRATOS_CREATE_3D_VARIABLE_WITH_THIS_COMPONENTS
#define KRATOS_CREATE_3D_VARIABLE_WITH_THIS_COMPONENTS(name, component1, component2, component3) ;

#undef KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS
#define KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(name) \
     KRATOS_CREATE_3D_VARIABLE_WITH_THIS_COMPONENTS(name, name##_X, name##_Y, name##_Z)

#undef KRATOS_CREATE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS
#define KRATOS_CREATE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, component1, component2, component3) ;

#undef KRATOS_CREATE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#define KRATOS_CREATE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
     KRATOS_CREATE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, name##_XX, name##_YY, name##_XY)

#undef KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS
#define KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, component1, component2, component3, component4, component5, component6) ;

#undef KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#define KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
     KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, name##_XX, name##_YY, name##_ZZ, name##_XY, name##_YZ, name##_XZ)

#undef KRATOS_CREATE_2D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS
#define KRATOS_CREATE_2D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, component1, component2, component3, component4) ;

#undef KRATOS_CREATE_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#define KRATOS_CREATE_2D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
     KRATOS_CREATE_2D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, name##_XX, name##_XY, name##_YX, name##_YY)

#undef KRATOS_CREATE_3D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS
#define KRATOS_CREATE_3D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, component1, component2, component3, component4, component5, component6, component7, component8, component9) ;

#undef KRATOS_CREATE_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#define KRATOS_CREATE_3D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
     KRATOS_CREATE_3D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, name##_XX, name##_XY, name##_XZ, name##_YX, name##_YY, name##_YZ, name##_ZX, name##_ZY, name##_ZZ)

#ifdef KRATOS_REGISTER_VARIABLE
#undef KRATOS_REGISTER_VARIABLE
#endif
#define KRATOS_REGISTER_VARIABLE(name) \
    AddKratosComponent(name.Name(), name); \
    KratosComponents<VariableData>::Add(name.Name(), name); \
    name.Register();

#ifdef KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_REGISTER_VARIABLE(name) \
    KRATOS_REGISTER_VARIABLE(name##_X) \
    KRATOS_REGISTER_VARIABLE(name##_Y) \
    KRATOS_REGISTER_VARIABLE(name##_Z)

#ifdef KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_REGISTER_VARIABLE(name) \
    KRATOS_REGISTER_VARIABLE(name##_XX) \
    KRATOS_REGISTER_VARIABLE(name##_YY) \
    KRATOS_REGISTER_VARIABLE(name##_XY)

#ifdef KRATOS_REGISTER_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_REGISTER_VARIABLE(name) \
    KRATOS_REGISTER_VARIABLE(name##_XX) \
    KRATOS_REGISTER_VARIABLE(name##_YY) \
    KRATOS_REGISTER_VARIABLE(name##_ZZ) \
    KRATOS_REGISTER_VARIABLE(name##_XY) \
    KRATOS_REGISTER_VARIABLE(name##_YZ) \
    KRATOS_REGISTER_VARIABLE(name##_XZ)

#ifdef KRATOS_REGISTER_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_2D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_REGISTER_VARIABLE(name) \
    KRATOS_REGISTER_VARIABLE(name##_XX) \
    KRATOS_REGISTER_VARIABLE(name##_XY) \
    KRATOS_REGISTER_VARIABLE(name##_YX) \
    KRATOS_REGISTER_VARIABLE(name##_YY)

#ifdef KRATOS_REGISTER_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_3D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_REGISTER_VARIABLE(name) \
    KRATOS_REGISTER_VARIABLE(name##_XX) \
    KRATOS_REGISTER_VARIABLE(name##_XY) \
    KRATOS_REGISTER_VARIABLE(name##_XZ) \
    KRATOS_REGISTER_VARIABLE(name##_YX) \
    KRATOS_REGISTER_VARIABLE(name##_YY) \
    KRATOS_REGISTER_VARIABLE(name##_YZ) \
    KRATOS_REGISTER_VARIABLE(name##_ZX) \
    KRATOS_REGISTER_VARIABLE(name##_ZY) \
    KRATOS_REGISTER_VARIABLE(name##_ZZ)

//-----------------------------------------------------------------
//
//  Variables time derivatives
//
//-----------------------------------------------------------------

#undef KRATOS_DEFINE_VARIABLE_WITH_TIME_DERIVATIVE
#define KRATOS_DEFINE_VARIABLE_WITH_TIME_DERIVATIVE(type, name, variable_derivative) \
    inline constexpr Variable<type> name = Variable<type>(#name, &variable_derivative);

#undef KRATOS_DEFINE_APPLICATION_VARIABLE_WITH_TIME_DERIVATIVE
#define KRATOS_DEFINE_APPLICATION_VARIABLE_WITH_TIME_DERIVATIVE(application, type, name, variable_derivative) \
    inline constexpr Variable<type> name = Variable<type>(#name, &variable_derivative);

#undef KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE
#define KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE(name, variable_derivative) \
    inline constexpr Variable<Kratos::array_1d<double, 3>> name = Variable<Kratos::array_1d<double, 3>>(#name, &variable_derivative); \
    inline constexpr Variable<double> CAT_(name, X) = Variable<double>(CAT_STR(name, X), &name, 0, &CAT_(variable_derivative, X)); \
    inline constexpr Variable<double> CAT_(name, Y) = Variable<double>(CAT_STR(name, Y), &name, 1, &CAT_(variable_derivative, Y)); \
    inline constexpr Variable<double> CAT_(name, Z) = Variable<double>(CAT_STR(name, Z), &name, 2, &CAT_(variable_derivative, Z));

#undef KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE
#define KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE(application, name, variable_derivative) \
    inline constexpr Variable<Kratos::array_1d<double, 3>> name = Variable<Kratos::array_1d<double, 3>>(#name, &variable_derivative); \
    inline constexpr Variable<double> CAT_(name, X) = Variable<double>(CAT_STR(name, X), &name, 0, &CAT_(variable_derivative, X)); \
    inline constexpr Variable<double> CAT_(name, Y) = Variable<double>(CAT_STR(name, Y), &name, 1, &CAT_(variable_derivative, Y)); \
    inline constexpr Variable<double> CAT_(name, Z) = Variable<double>(CAT_STR(name, Z), &name, 2, &CAT_(variable_derivative, Z));

#undef KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE
#define KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE(name, variable_derivative) \
    inline constexpr Variable<Kratos::array_1d<double, 3>> name = Variable<Kratos::array_1d<double, 3>>(#name, &variable_derivative); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0, &CAT_(variable_derivative, XX)); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, YY), &name, 1, &CAT_(variable_derivative, YY)); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 2, &CAT_(variable_derivative, XY));

#undef KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE
#define KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE(application, name, variable_derivative) \
    inline constexpr Variable<Kratos::array_1d<double, 3>> name = Variable<Kratos::array_1d<double, 3>>(#name, &variable_derivative); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0, &CAT_(variable_derivative, XX)); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, YY), &name, 1, &CAT_(variable_derivative, YY)); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 2, &CAT_(variable_derivative, XY));

#undef KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE
#define KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE(name, variable_derivative) \
    inline constexpr Variable<Kratos::array_1d<double, 6>> name = Variable<Kratos::array_1d<double, 6>>(#name, &variable_derivative); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0, &CAT_(variable_derivative, XX)); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, YY), &name, 1, &CAT_(variable_derivative, YY)); \
    inline constexpr Variable<double> CAT_(name, ZZ) = Variable<double>(CAT_STR(name, ZZ), &name, 2, &CAT_(variable_derivative, ZZ)); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 3, &CAT_(variable_derivative, XY)); \
    inline constexpr Variable<double> CAT_(name, YZ) = Variable<double>(CAT_STR(name, YZ), &name, 4, &CAT_(variable_derivative, YZ)); \
    inline constexpr Variable<double> CAT_(name, XZ) = Variable<double>(CAT_STR(name, XZ), &name, 5, &CAT_(variable_derivative, XZ));

#undef KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE
#define KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE(application, name, variable_derivative) \
    inline constexpr Variable<Kratos::array_1d<double, 6>> name = Variable<Kratos::array_1d<double, 6>>(#name, &variable_derivative); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0, &CAT_(variable_derivative, XX)); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, YY), &name, 1, &CAT_(variable_derivative, YY)); \
    inline constexpr Variable<double> CAT_(name, ZZ) = Variable<double>(CAT_STR(name, ZZ), &name, 2, &CAT_(variable_derivative, ZZ)); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 3, &CAT_(variable_derivative, XY)); \
    inline constexpr Variable<double> CAT_(name, YZ) = Variable<double>(CAT_STR(name, YZ), &name, 4, &CAT_(variable_derivative, YZ)); \
    inline constexpr Variable<double> CAT_(name, XZ) = Variable<double>(CAT_STR(name, XZ), &name, 5, &CAT_(variable_derivative, XZ));

#undef KRATOS_DEFINE_2D_TENSOR_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE
#define KRATOS_DEFINE_2D_TENSOR_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE(name, variable_derivative) \
    inline constexpr Variable<Kratos::array_1d<double, 4>> name = Variable<Kratos::array_1d<double, 4>>(#name, &variable_derivative); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0, &CAT_(variable_derivative, XX)); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 1, &CAT_(variable_derivative, XY)); \
    inline constexpr Variable<double> CAT_(name, YX) = Variable<double>(CAT_STR(name, YX), &name, 2, &CAT_(variable_derivative, YX)); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, YY), &name, 3, &CAT_(variable_derivative, YY));

#undef KRATOS_DEFINE_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE
#define KRATOS_DEFINE_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE(application, name, variable_derivative) \
    inline constexpr Variable<Kratos::array_1d<double, 4>> name = Variable<Kratos::array_1d<double, 4>>(#name, &variable_derivative); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0, &CAT_(variable_derivative, XX)); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 1, &CAT_(variable_derivative, XY)); \
    inline constexpr Variable<double> CAT_(name, YX) = Variable<double>(CAT_STR(name, YX), &name, 2, &CAT_(variable_derivative, YX)); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, YY), &name, 3, &CAT_(variable_derivative, YY));

#undef KRATOS_DEFINE_3D_TENSOR_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE
#define KRATOS_DEFINE_3D_TENSOR_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE(name, variable_derivative) \
    inline constexpr Variable<Kratos::array_1d<double, 9>> name = Variable<Kratos::array_1d<double, 9>>(#name, &variable_derivative); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0, &CAT_(variable_derivative, XX)); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 1, &CAT_(variable_derivative, XY)); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 2, &CAT_(variable_derivative, XY)); \
    inline constexpr Variable<double> CAT_(name, YX) = Variable<double>(CAT_STR(name, YX), &name, 3, &CAT_(variable_derivative, YX)); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, YY), &name, 4, &CAT_(variable_derivative, YY)); \
    inline constexpr Variable<double> CAT_(name, YZ) = Variable<double>(CAT_STR(name, YZ), &name, 5, &CAT_(variable_derivative, YZ)); \
    inline constexpr Variable<double> CAT_(name, ZX) = Variable<double>(CAT_STR(name, ZX), &name, 6, &CAT_(variable_derivative, ZX)); \
    inline constexpr Variable<double> CAT_(name, ZY) = Variable<double>(CAT_STR(name, ZY), &name, 7, &CAT_(variable_derivative, ZY)); \
    inline constexpr Variable<double> CAT_(name, ZZ) = Variable<double>(CAT_STR(name, ZZ), &name, 8, &CAT_(variable_derivative, ZZ)); \

#undef KRATOS_DEFINE_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE
#define KRATOS_DEFINE_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS_WITH_TIME_DERIVATIVE(application, name, variable_derivative) \
    inline constexpr Variable<Kratos::array_1d<double, 9>> name = Variable<Kratos::array_1d<double, 9>>(#name, &variable_derivative); \
    inline constexpr Variable<double> CAT_(name, XX) = Variable<double>(CAT_STR(name, XX), &name, 0, &CAT_(variable_derivative, XX)); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 1, &CAT_(variable_derivative, XY)); \
    inline constexpr Variable<double> CAT_(name, XY) = Variable<double>(CAT_STR(name, XY), &name, 2, &CAT_(variable_derivative, XY)); \
    inline constexpr Variable<double> CAT_(name, YX) = Variable<double>(CAT_STR(name, YX), &name, 3, &CAT_(variable_derivative, YX)); \
    inline constexpr Variable<double> CAT_(name, YY) = Variable<double>(CAT_STR(name, YY), &name, 4, &CAT_(variable_derivative, YY)); \
    inline constexpr Variable<double> CAT_(name, YZ) = Variable<double>(CAT_STR(name, YZ), &name, 5, &CAT_(variable_derivative, YZ)); \
    inline constexpr Variable<double> CAT_(name, ZX) = Variable<double>(CAT_STR(name, ZX), &name, 6, &CAT_(variable_derivative, ZX)); \
    inline constexpr Variable<double> CAT_(name, ZY) = Variable<double>(CAT_STR(name, ZY), &name, 7, &CAT_(variable_derivative, ZY)); \
    inline constexpr Variable<double> CAT_(name, ZZ) = Variable<double>(CAT_STR(name, ZZ), &name, 8, &CAT_(variable_derivative, ZZ)); \

//-----------------------------------------------------------------
//
// Flags
//
//-----------------------------------------------------------------

#ifdef KRATOS_DEFINE_FLAG
#undef KRATOS_DEFINE_FLAG
#endif
#define KRATOS_DEFINE_FLAG(name) \
    extern const Kratos::Flags name;

#ifdef KRATOS_ADD_FLAG_TO_KRATOS_COMPONENTS
#undef KRATOS_ADD_FLAG_TO_KRATOS_COMPONENTS
#endif
#define KRATOS_ADD_FLAG_TO_KRATOS_COMPONENTS(name)                  \
    Kratos::KratosComponents<Kratos::Flags>::Add(#name, name)

#ifdef KRATOS_CREATE_FLAG
#undef KRATOS_CREATE_FLAG
#endif
#define KRATOS_CREATE_FLAG(name, position)                  \
    const Kratos::Flags name(Kratos::Flags::Create(position));

#ifdef KRATOS_REGISTER_FLAG
#undef KRATOS_REGISTER_FLAG
#endif
#define KRATOS_REGISTER_FLAG(name)                  \
    KRATOS_ADD_FLAG_TO_KRATOS_COMPONENTS(name);



#ifdef KRATOS_DEFINE_LOCAL_FLAG
#undef KRATOS_DEFINE_LOCAL_FLAG
#endif
#define KRATOS_DEFINE_LOCAL_FLAG(name)		\
  static const Kratos::Flags name;

#ifdef KRATOS_DEFINE_LOCAL_APPLICATION_FLAG
#undef KRATOS_DEFINE_LOCAL_APPLICATION_FLAG
#endif
#define KRATOS_DEFINE_LOCAL_APPLICATION_FLAG(application, name)		\
  static const Kratos::Flags name;

#ifdef KRATOS_CREATE_LOCAL_FLAG
#undef KRATOS_CREATE_LOCAL_FLAG
#endif
#define KRATOS_CREATE_LOCAL_FLAG(class_name, name, position)		\
  const Kratos::Flags class_name::name(Kratos::Flags::Create(position));



//-----------------------------------------------------------------
//
// components
//
//-----------------------------------------------------------------

#ifdef KRATOS_REGISTER_GEOMETRY
#undef KRATOS_REGISTER_GEOMETRY
#endif
#define KRATOS_REGISTER_GEOMETRY(name, reference)                                                                   \
    KratosComponents<Geometry<Node>>::Add(name, reference);                                                         \
    if(!Registry::HasItem("geometries."+Registry::GetCurrentSource()+"."+name) &&                                   \
       !Registry::HasItem("components."+std::string(name))){                                                                     \
        Registry::AddItem<RegistryItem>("geometries."+Registry::GetCurrentSource()+"."+name);                       \
        Registry::AddItem<RegistryItem>("components."+std::string(name));                                                        \
    }                                                                                                               \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_ELEMENT
#undef KRATOS_REGISTER_ELEMENT
#endif
#define KRATOS_REGISTER_ELEMENT(name, reference)                                                                    \
    KratosComponents<Element>::Add(name, reference);                                                                \
    if(!Registry::HasItem("elements."+Registry::GetCurrentSource()+"."+name) &&                                     \
       !Registry::HasItem("components."+std::string(name))){                                                                    \
        Registry::AddItem<RegistryItem>("elements."+Registry::GetCurrentSource()+"."+name);                         \
        Registry::AddItem<RegistryItem>("components."+std::string(name));                                                        \
    }                                                                                                               \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_CONDITION
#undef KRATOS_REGISTER_CONDITION
#endif
#define KRATOS_REGISTER_CONDITION(name, reference)                                                                  \
    KratosComponents<Condition>::Add(name, reference);                                                              \
    if(!Registry::HasItem("conditions."+Registry::GetCurrentSource()+"."+name) &&                                   \
       !Registry::HasItem("components."+std::string(name))){                                                                     \
        Registry::AddItem<RegistryItem>("conditions."+Registry::GetCurrentSource()+"."+name);                       \
        Registry::AddItem<RegistryItem>("components."+std::string(name));                                                        \
    }                                                                                                               \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_CONSTRAINT
#undef KRATOS_REGISTER_CONSTRAINT
#endif
#define KRATOS_REGISTER_CONSTRAINT(name, reference)                                                                 \
    KratosComponents<MasterSlaveConstraint>::Add(name, reference);                                                  \
    if(!Registry::HasItem("constraints."+Registry::GetCurrentSource()+"."+name) &&                                  \
       !Registry::HasItem("components."+std::string(name))){                                                                     \
        Registry::AddItem<RegistryItem>("constraints."+Registry::GetCurrentSource()+"."+name);                      \
        Registry::AddItem<RegistryItem>("components."+std::string(name));                                                        \
    }                                                                                                               \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_MODELER
#undef KRATOS_REGISTER_MODELER
#endif
#define KRATOS_REGISTER_MODELER(name, reference)                                                                    \
    KratosComponents<Modeler>::Add(name, reference);                                                                \
    if(!Registry::HasItem("modelers."+Registry::GetCurrentSource()+"."+name) &&                                     \
       !Registry::HasItem("components."+std::string(name))){                                                                     \
        Registry::AddItem<RegistryItem>("modelers."+Registry::GetCurrentSource()+"."+name);                         \
        Registry::AddItem<RegistryItem>("components."+std::string(name));                                                        \
    }                                                                                                               \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_CONSTITUTIVE_LAW
#undef KRATOS_REGISTER_CONSTITUTIVE_LAW
#endif
#define KRATOS_REGISTER_CONSTITUTIVE_LAW(name, reference)                                                           \
    KratosComponents<ConstitutiveLaw>::Add(name, reference);                                                        \
    if(!Registry::HasItem("constitutive_laws."+Registry::GetCurrentSource()+"."+name) &&                            \
       !Registry::HasItem("components."+std::string(name))){                                                                     \
        Registry::AddItem<RegistryItem>("constitutive_laws."+Registry::GetCurrentSource()+"."+name);                \
        Registry::AddItem<RegistryItem>("components."+std::string(name));                                                        \
    }                                                                                                               \
    Serializer::Register(name, reference);

#define KRATOS_DEPRECATED [[deprecated]]
#define KRATOS_DEPRECATED_MESSAGE(deprecated_message) [[deprecated(deprecated_message)]]

// The following block defines the macro KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
// If written in a file, for the following lines of code the compiler will not print warnings of type 'deprecated function'.
// The scope ends where KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING is called.
// NOTE!! this macro is not intended for extensive use, it's just for temporary use in methods exported to Python which
// are still calling a C++ deprecated function.
#if defined(__clang__)
#define KRATOS_PRAGMA_INSIDE_MACRO_DEFINITION(x) _Pragma(#x)
#define KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING \
KRATOS_PRAGMA_INSIDE_MACRO_DEFINITION(clang diagnostic push) \
KRATOS_PRAGMA_INSIDE_MACRO_DEFINITION(clang diagnostic ignored "-Wdeprecated-declarations")
#elif defined(__GNUG__) && !defined(__INTEL_COMPILER)
#define KRATOS_PRAGMA_INSIDE_MACRO_DEFINITION(x) _Pragma(#x)
#define KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING \
KRATOS_PRAGMA_INSIDE_MACRO_DEFINITION(GCC diagnostic push) \
KRATOS_PRAGMA_INSIDE_MACRO_DEFINITION(GCC diagnostic ignored "-Wdeprecated-declarations")
#elif defined(_MSC_VER)
#define KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING \
__pragma(warning(push))\
__pragma(warning(disable: 4996))
#else
#define KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING // not implemented for other compilers, hence does nothing
#endif

// The following block defines the macro KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING which ends the scope for
// ignoring the warnings of type 'deprecated function'.
#if defined(__clang__)
#define KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING \
_Pragma("clang diagnostic pop")
#elif defined(__GNUG__) && !defined(__INTEL_COMPILER)
#define KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING \
_Pragma("GCC diagnostic pop")
#elif defined(_MSC_VER)
#define KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING \
__pragma(warning(pop))
#else
#define KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING // not implemented for other compilers, hence does nothing
#endif


namespace Kratos
{
///@name Type Definitions
///@{

#if defined(_MSC_VER)
#pragma warning (disable: 4355)
#pragma warning (disable: 4503)
#pragma warning (disable: 4786)
#endif

//Exception handling
#define KRATOS_TYPE_NAME_OF(name) name##Type
#define KRATOS_NOT_EXCLUDED(filename) !defined(KRATOS_##filename##_EXCLUDED)

#define KRATOS_DECLEAR_TYPE  namespace KratosComponents{ typedef
#define KRATOS_FOR_COMPONENT_NAMED(name) KRATOS_TYPE_NAME_OF(name);}

// Kratos variable registering
/* #define KRATOS_REGISTER_VARIABLE_WITH_ZERO(type, name, zero) const Variable<type > name(#name, __LINE__, zero) */
/* #define KRATOS_REGISTER_VARIABLE(type, name) const Variable<type > name(#name, __LINE__) */

/* #define KRATOS_REGISTER_LINEAR_SOLVER_BEGIN \ */
/* template<class TFunction> ApplyToLinearSolver(String Name){ */

//Print Trace if defined
#define KRATOS_WATCH(variable) std::cout << #variable << " : " << variable << std::endl;
#define KRATOS_WATCH_CERR(variable) std::cerr << #variable << " : " << variable << std::endl;
#define KRATOS_WATCH_MPI(variable, mpi_data_comm) std::cout << "RANK " << mpi_data_comm.Rank() << "/" << mpi_data_comm.Size()  << "    "; KRATOS_WATCH(variable);

}  /* namespace Kratos.*/

#define KRATOS_SERIALIZE_SAVE_BASE_CLASS(Serializer, BaseType) \
    Serializer.save_base("BaseClass",*static_cast<const BaseType *>(this));

#define KRATOS_SERIALIZE_LOAD_BASE_CLASS(Serializer, BaseType) \
    Serializer.load_base("BaseClass",*static_cast<BaseType *>(this));
