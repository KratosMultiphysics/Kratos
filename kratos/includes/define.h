//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#if !defined(KRATOS_DEFINE_H_INCLUDED )
#define  KRATOS_DEFINE_H_INCLUDED

/* System includes */
#include <stdexcept>
#include <sstream>


/* External includes */

/* Project includes */
#include "includes/kratos_export_api.h"
#include "includes/shared_pointers.h"
#include "includes/exception.h"


#if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    #define KRATOS_COMPILED_IN_LINUX

#elif defined(__APPLE__) && defined(__MACH__)
    #define KRATOS_COMPILED_IN_OSX

#elif defined(_WIN32)
    #define KRATOS_COMPILED_IN_WINDOWS
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

#ifndef __SUNPRO_CC
#define KRATOS_TRY try {

#define KRATOS_CATCH(MoreInfo) \
  KRATOS_CATCH_WITH_BLOCK(MoreInfo,{})
#else
#define KRATOS_TRY { };

#define KRATOS_CATCH(MoreInfo) { };
#endif

//-----------------------------------------------------------------
//
// variables
//
//-----------------------------------------------------------------

#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

#ifdef KRATOS_DEFINE_VARIABLE_IMPLEMENTATION
#undef KRATOS_DEFINE_VARIABLE_IMPLEMENTATION
#endif
#define KRATOS_DEFINE_VARIABLE_IMPLEMENTATION(module, type, name) \
    KRATOS_EXPORT_MACRO(module) extern Variable<type > name;

#ifdef KRATOS_DEFINE_VARIABLE
#undef KRATOS_DEFINE_VARIABLE
#endif
#define KRATOS_DEFINE_VARIABLE(type, name) \
    KRATOS_DEFINE_VARIABLE_IMPLEMENTATION(KRATOS_CORE, type, name)

#ifdef KRATOS_DEFINE_APPLICATION_VARIABLE
#undef KRATOS_DEFINE_APPLICATION_VARIABLE
#endif
#define KRATOS_DEFINE_APPLICATION_VARIABLE(application, type, name) \
    KRATOS_API(application) extern Variable<type > name;

#ifdef KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS_IMPLEMENTATION(module, name) \
    KRATOS_EXPORT_MACRO(module) extern Kratos::Variable<Kratos::array_1d<double, 3> > name; \
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_X;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_Y;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_Z;

#ifdef KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS_IMPLEMENTATION(KRATOS_CORE, name)

#ifdef KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(application, name) \
  KRATOS_API(application) extern Kratos::Variable<Kratos::array_1d<double, 3> > name; \
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_X;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_Y;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_Z;

#ifdef KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS_IMPLEMENTATION(module, name) \
    KRATOS_EXPORT_MACRO(module) extern Kratos::Variable<Kratos::array_1d<double, 3> > name; \
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_XX;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_YY;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_XY;

#ifdef KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS_IMPLEMENTATION(KRATOS_CORE, name)

#ifdef KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS(application, name) \
  KRATOS_API(application) extern Kratos::Variable<Kratos::array_1d<double, 3> > name; \
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_XX;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_YY;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_XY;

#ifdef KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS_IMPLEMENTATION(module, name) \
    KRATOS_EXPORT_MACRO(module) extern Kratos::Variable<Kratos::array_1d<double, 6> > name; \
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > name##_XX;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > name##_YY;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > name##_ZZ;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > name##_XY;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > name##_YZ;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > name##_XZ;

#ifdef KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS_IMPLEMENTATION(KRATOS_CORE, name)

#ifdef KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS(application, name) \
  KRATOS_API(application) extern Kratos::Variable<Kratos::array_1d<double, 6> > name; \
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > name##_XX;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > name##_YY;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > name##_ZZ;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > name##_XY;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > name##_YZ;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > name##_XZ;

#ifdef KRATOS_DEFINE_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_2D_TENSOR_VARIABLE_WITH_COMPONENTS_IMPLEMENTATION(module, name) \
    KRATOS_EXPORT_MACRO(module) extern Kratos::Variable<Kratos::array_1d<double, 4> > name; \
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> > > name##_XX;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> > > name##_XY;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> > > name##_YX;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> > > name##_YY;

#ifdef KRATOS_DEFINE_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_2D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_DEFINE_2D_TENSOR_VARIABLE_WITH_COMPONENTS_IMPLEMENTATION(KRATOS_CORE, name)

#ifdef KRATOS_DEFINE_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS(application, name) \
  KRATOS_API(application) extern Kratos::Variable<Kratos::array_1d<double, 4> > name; \
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> > > name##_XX;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> > > name##_XY;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> > > name##_YX;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> > > name##_YY;

#ifdef KRATOS_DEFINE_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_3D_TENSOR_VARIABLE_WITH_COMPONENTS_IMPLEMENTATION(module, name) \
    KRATOS_EXPORT_MACRO(module) extern Kratos::Variable<Kratos::array_1d<double, 9> > name; \
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_XX;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_XY;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_XZ;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_YX;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_YY;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_YZ;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_ZX;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_ZY;\
    KRATOS_EXPORT_MACRO(module) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_ZZ;

#ifdef KRATOS_DEFINE_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_3D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_DEFINE_3D_TENSOR_VARIABLE_WITH_COMPONENTS_IMPLEMENTATION(KRATOS_CORE, name)

#ifdef KRATOS_DEFINE_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS(application, name) \
  KRATOS_API(application) extern Kratos::Variable<Kratos::array_1d<double, 9> > name; \
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_XX;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_XY;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_XZ;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_YX;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_YY;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_YZ;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_ZX;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_ZY;\
  KRATOS_API(application) extern Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > name##_ZZ;

#ifdef KRATOS_CREATE_VARIABLE
#undef KRATOS_CREATE_VARIABLE
#endif
#define KRATOS_CREATE_VARIABLE(type, name) \
    /*const*/ Kratos::Variable<type > name(#name);

#ifdef KRATOS_CREATE_VARIABLE_WITH_ZERO
#undef KRATOS_CREATE_VARIABLE_WITH_ZERO
#endif
#define KRATOS_CREATE_VARIABLE_WITH_ZERO(type, name, zero) \
    /*const*/ Kratos::Variable<type> name(#name, zero);

#ifdef KRATOS_CREATE_3D_VARIABLE_WITH_THIS_COMPONENTS
#undef KRATOS_CREATE_3D_VARIABLE_WITH_THIS_COMPONENTS
#endif
#define KRATOS_CREATE_3D_VARIABLE_WITH_THIS_COMPONENTS(name, component1, component2, component3) \
    /*const*/ Kratos::Variable<Kratos::array_1d<double, 3> > name(#name, Kratos::array_1d<double, 3>(Kratos::ZeroVector(3))); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > \
                  component1(#component1, #name, 0, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >(name, 0)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > \
                  component2(#component2, #name, 1, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >(name, 1)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > \
                  component3(#component3, #name, 2, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >(name, 2));

#ifdef KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS
#undef KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(name) \
     KRATOS_CREATE_3D_VARIABLE_WITH_THIS_COMPONENTS(name, name##_X, name##_Y, name##_Z)

#ifdef KRATOS_CREATE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS
#undef KRATOS_CREATE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS
#endif
#define KRATOS_CREATE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, component1, component2, component3) \
    /*const*/ Kratos::Variable<Kratos::array_1d<double, 3> > name(#name, Kratos::zero_vector<double>(3)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > \
                  component1(#component1, #name, 0, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >(name, 0)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > \
                  component2(#component2, #name, 1, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >(name, 1)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > \
                  component3(#component3, #name, 2, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >(name, 2));

#ifdef KRATOS_CREATE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_CREATE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_CREATE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
     KRATOS_CREATE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, name##_XX, name##_YY, name##_XY)

#ifdef KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS
#undef KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS
#endif
#define KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, component1, component2, component3, component4, component5, component6) \
    /*const*/ Kratos::Variable<Kratos::array_1d<double, 6> > name(#name, Kratos::zero_vector<double>(6)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > \
                  component1(#component1, #name, 0, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> >(name, 0)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > \
                  component2(#component2, #name, 1, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> >(name, 1)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > \
                  component3(#component3, #name, 2, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> >(name, 2)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > \
                  component4(#component4, #name, 3, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> >(name, 3)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > \
                  component5(#component5, #name, 4, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> >(name, 4)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> > > \
                  component6(#component6, #name, 5, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> >(name, 5));

#ifdef KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
     KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, name##_XX, name##_YY, name##_ZZ, name##_XY, name##_YZ, name##_XZ)

#ifdef KRATOS_CREATE_2D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS
#undef KRATOS_CREATE_2D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS
#endif
#define KRATOS_CREATE_2D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, component1, component2, component3, component4) \
    /*const*/ Kratos::Variable<Kratos::array_1d<double, 4> > name(#name, Kratos::zero_vector<double>(4)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> > > \
                  component1(#component1, #name, 0, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> >(name, 0)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> > > \
                  component2(#component2, #name, 1, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> >(name, 1)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> > > \
                  component3(#component3, #name, 2, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> >(name, 2)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> > > \
                  component4(#component4, #name, 3, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> >(name, 3));

#ifdef KRATOS_CREATE_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_CREATE_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_CREATE_2D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
     KRATOS_CREATE_2D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, name##_XX, name##_XY, name##_YX, name##_YY)

#ifdef KRATOS_CREATE_3D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS
#undef KRATOS_CREATE_3D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS
#endif
#define KRATOS_CREATE_3D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, component1, component2, component3, component4, component5, component6, component7, component8, component9) \
    /*const*/ Kratos::Variable<Kratos::array_1d<double, 9> > name(#name, Kratos::zero_vector<double>(9)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > \
                  component1(#component1, #name, 0, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> >(name, 0)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > \
                  component2(#component2, #name, 1, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> >(name, 1)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > \
                  component3(#component3, #name, 2, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> >(name, 2)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > \
                  component4(#component4, #name, 3, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> >(name, 3)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > \
                  component5(#component5, #name, 4, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> >(name, 4)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > \
                  component6(#component6, #name, 5, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> >(name, 5)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > \
                  component7(#component7, #name, 6, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> >(name, 6)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > \
                  component8(#component8, #name, 7, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> >(name, 7)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> > > \
                  component9(#component9, #name, 8, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> >(name, 8));

#ifdef KRATOS_CREATE_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_CREATE_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_CREATE_3D_TENSOR_VARIABLE_WITH_COMPONENTS(name) \
     KRATOS_CREATE_3D_TENSOR_VARIABLE_WITH_THIS_COMPONENTS(name, name##_XX, name##_XY, name##_XZ, name##_YX, name##_YY, name##_YZ, name##_ZX, name##_ZY, name##_ZZ)

#ifdef KRATOS_REGISTER_VARIABLE
#undef KRATOS_REGISTER_VARIABLE
#endif
#define KRATOS_REGISTER_VARIABLE(name) \
    AddKratosComponent(name.Name(), name); \
    KratosComponents<VariableData>::Add(name.Name(), name);

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
// Flags
//
//-----------------------------------------------------------------

#ifdef KRATOS_DEFINE_FLAG
#undef KRATOS_DEFINE_FLAG
#endif
#define KRATOS_DEFINE_FLAG(name) \
    extern const Kratos::Flags name;     \
    extern const Kratos::Flags NOT_##name

#ifdef KRATOS_ADD_FLAG_TO_KRATOS_COMPONENTS
#undef KRATOS_ADD_FLAG_TO_KRATOS_COMPONENTS
#endif
#define KRATOS_ADD_FLAG_TO_KRATOS_COMPONENTS(name)                  \
    Kratos::KratosComponents<Kratos::Flags>::Add(#name, name)

#ifdef KRATOS_CREATE_FLAG
#undef KRATOS_CREATE_FLAG
#endif
#define KRATOS_CREATE_FLAG(name, position)                  \
    const Kratos::Flags name(Kratos::Flags::Create(position));              \
    const Kratos::Flags NOT_##name(Kratos::Flags::Create(position, false))

#ifdef KRATOS_REGISTER_FLAG
#undef KRATOS_REGISTER_FLAG
#endif
#define KRATOS_REGISTER_FLAG(name)                  \
    KRATOS_ADD_FLAG_TO_KRATOS_COMPONENTS(name);             \
    KRATOS_ADD_FLAG_TO_KRATOS_COMPONENTS(NOT_##name)



#ifdef KRATOS_DEFINE_LOCAL_FLAG
#undef KRATOS_DEFINE_LOCAL_FLAG
#endif
#define KRATOS_DEFINE_LOCAL_FLAG(name)		\
  static const Kratos::Flags name;			\
  static const Kratos::Flags NOT_##name

#ifdef KRATOS_DEFINE_LOCAL_APPLICATION_FLAG
#undef KRATOS_DEFINE_LOCAL_APPLICATION_FLAG
#endif
#define KRATOS_DEFINE_LOCAL_APPLICATION_FLAG(application, name)		\
  static KRATOS_API(application) const Kratos::Flags name;			\
  static KRATOS_API(application) const Kratos::Flags NOT_##name

#ifdef KRATOS_CREATE_LOCAL_FLAG
#undef KRATOS_CREATE_LOCAL_FLAG
#endif
#define KRATOS_CREATE_LOCAL_FLAG(class_name, name, position)		\
  const Kratos::Flags class_name::name(Kratos::Flags::Create(position));		\
  const Kratos::Flags class_name::NOT_##name(Kratos::Flags::Create(position, false))



//-----------------------------------------------------------------
//
// components
//
//-----------------------------------------------------------------

#ifdef KRATOS_REGISTER_ELEMENT
#undef KRATOS_REGISTER_ELEMENT
#endif
#define KRATOS_REGISTER_ELEMENT(name, reference) \
    KratosComponents<Element >::Add(name, reference); \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_CONDITION
#undef KRATOS_REGISTER_CONDITION
#endif
#define KRATOS_REGISTER_CONDITION(name, reference) \
    KratosComponents<Condition >::Add(name, reference); \
    Serializer::Register(name, reference);

#ifdef KRATOS_REGISTER_CONSTRAINT
#undef KRATOS_REGISTER_CONSTRAINT
#endif
#define KRATOS_REGISTER_CONSTRAINT(name, reference) \
    KratosComponents<MasterSlaveConstraint >::Add(name, reference); \
    Serializer::Register(name, reference);


#ifdef KRATOS_REGISTER_CONSTITUTIVE_LAW
#undef KRATOS_REGISTER_CONSTITUTIVE_LAW
#endif
#define KRATOS_REGISTER_CONSTITUTIVE_LAW(name, reference) \
    KratosComponents<ConstitutiveLaw >::Add(name, reference); \
    Serializer::Register(name, reference);

#if __cplusplus >= 201402L
#define KRATOS_DEPRECATED [[deprecated]]
#define KRATOS_DEPRECATED_MESSAGE(deprecated_message) [[deprecated(deprecated_message)]]
#elif __GNUC__
#define KRATOS_DEPRECATED __attribute__((deprecated))
#define KRATOS_DEPRECATED_MESSAGE(deprecated_message) KRATOS_DEPRECATED
#elif defined(_MSC_VER)
#define KRATOS_DEPRECATED __declspec(deprecated)
#define KRATOS_DEPRECATED_MESSAGE(deprecated_message) KRATOS_DEPRECATED
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define KRATOS_DEPRECATED
#define KRATOS_DEPRECATED_MESSAGE(deprecated_message)
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
/* #define KRATOS_REGISTER_VARIABLE_COMPONENT(type, name, source) const VariableComponent<type > name(#name, __LINE__, type source) */

/* #define KRATOS_REGISTER_LINEAR_SOLVER_BEGIN \ */
/* template<class TFunction> ApplyToLinearSolver(String Name){ */

//Print Trace if defined
#define KRATOS_WATCH(variable) std::cout << #variable << " : " << variable << std::endl;

}  /* namespace Kratos.*/

#define KRATOS_SERIALIZE_SAVE_BASE_CLASS(Serializer, BaseType) \
	Serializer.save_base("BaseClass",*static_cast<const BaseType *>(this));

#define KRATOS_SERIALIZE_LOAD_BASE_CLASS(Serializer, BaseType) \
	Serializer.load_base("BaseClass",*static_cast<BaseType *>(this));


#endif /* KRATOS_DEFINE_H_INCLUDED  defined */
