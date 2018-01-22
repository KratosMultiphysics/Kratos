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
#include "boost/smart_ptr.hpp"
#include "boost/current_function.hpp"


/* Project includes */
#include "includes/kratos_config.h"
#include "includes/kratos_export_api.h"
#include "includes/exception.h"



namespace Kratos
{
template<class T>
    using shared_ptr = boost::shared_ptr<T>; //std::shared_ptr<T>;

template<class T>
    using weak_ptr = boost::weak_ptr<T>; //std::weak_ptr<T>;
    
template<typename C, typename...Args>
    shared_ptr<C> make_shared(Args &&...args) {
        return boost::make_shared<C>(std::forward<Args>(args)...);
//         return std::make_shared<C>(std::forward<Args>(args)...);
    }

}
    


#define KRATOS_CLASS_POINTER_DEFINITION(a) typedef Kratos::shared_ptr<a > Pointer; \
typedef Kratos::shared_ptr<a > SharedPointer; \
typedef Kratos::weak_ptr<a > WeakPointer


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

#if defined(KRATOS_SET_EXCEPTION_LEVEL_TO_1)
#define KRATOS_EXCEPTION_LEVEL_1
#endif

#if defined(KRATOS_SET_EXCEPTION_LEVEL_TO_2)
#define KRATOS_EXCEPTION_LEVEL_1
#define KRATOS_EXCEPTION_LEVEL_2
#endif

#if defined(KRATOS_SET_EXCEPTION_LEVEL_TO_3)
#define KRATOS_EXCEPTION_LEVEL_1
#define KRATOS_EXCEPTION_LEVEL_2
#define KRATOS_EXCEPTION_LEVEL_3
#endif

#if defined(KRATOS_SET_EXCEPTION_LEVEL_TO_4)
#define KRATOS_EXCEPTION_LEVEL_1
#define KRATOS_EXCEPTION_LEVEL_2
#define KRATOS_EXCEPTION_LEVEL_3
#define KRATOS_EXCEPTION_LEVEL_4
#endif

#if defined(KRATOS_EXCEPTION_LEVEL_1)
#define KRATOS_TRY_LEVEL_1 try {
#define KRATOS_CATCH_LEVEL_1(MoreInfo) \
KRATOS_CATCH_WITH_BLOCK(MoreInfo,{})
#else
#define KRATOS_TRY_LEVEL_1 {
#define KRATOS_CATCH_LEVEL_1(MoreInfo) }
#endif


#if defined(KRATOS_EXCEPTION_LEVEL_2)
#define KRATOS_TRY_LEVEL_2 try {
#define KRATOS_CATCH_LEVEL_2(MoreInfo) \
KRATOS_CATCH_WITH_BLOCK(MoreInfo,{})
#else
#define KRATOS_TRY_LEVEL_2 {
#define KRATOS_CATCH_LEVEL_2(MoreInfo) }
#endif

#if defined(KRATOS_EXCEPTION_LEVEL_3)
#define KRATOS_TRY_LEVEL_3 try {
#define KRATOS_CATCH_LEVEL_3(MoreInfo) \
KRATOS_CATCH_WITH_BLOCK(MoreInfo,{})
#else
#define KRATOS_TRY_LEVEL_3 {
#define KRATOS_CATCH_LEVEL_3(MoreInfo) }
#endif

#if defined(KRATOS_EXCEPTION_LEVEL_4)
#define KRATOS_TRY_LEVEL_4 try {
#define KRATOS_CATCH_LEVEL_4(MoreInfo) \
KRATOS_CATCH_WITH_BLOCK(MoreInfo,{})
#else
#define KRATOS_TRY_LEVEL_4 {
#define KRATOS_CATCH_LEVEL_4(MoreInfo) }
#endif

#ifndef KRATOS_CURRENT_FUNCTION
#define KRATOS_CURRENT_FUNCTION BOOST_CURRENT_FUNCTION
#endif


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
    /*const*/ Kratos::Variable<Kratos::array_1d<double, 3> > name(#name, Kratos::zero_vector<double>(3)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > \
                  component1(#component1, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >(name, 0)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > \
                  component2(#component2, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >(name, 1)); \
\
    /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > \
                  component3(#component3, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >(name, 2));

#ifdef KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS
#undef KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(name) \
     KRATOS_CREATE_3D_VARIABLE_WITH_THIS_COMPONENTS(name, name##_X, name##_Y, name##_Z)

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

#ifdef KRATOS_REGISTER_IN_PYTHON_VARIABLE
#undef KRATOS_REGISTER_IN_PYTHON_VARIABLE
#endif
#define KRATOS_REGISTER_IN_PYTHON_VARIABLE(variable) \
	scope().attr(#variable) = boost::ref(variable);

#ifdef KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name##_X) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name##_Y) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name##_Z)

#ifdef KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION
#undef KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION
#endif
#define KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION(flag) \
    scope().attr(#flag) = boost::ref(flag)      \
 
#ifdef KRATOS_REGISTER_IN_PYTHON_FLAG
#undef KRATOS_REGISTER_IN_PYTHON_FLAG
#endif
#define KRATOS_REGISTER_IN_PYTHON_FLAG(flag) \
    KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION(flag);   \
    KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION(NOT_##flag)

    
    
    
#ifdef __GNUC__
#define KRATOS_DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define KRATOS_DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define KRATOS_DEPRECATED
#endif
    

namespace Kratos
{


/**@name Kratos Classes */
/*@{ */

/*@} */

/**@name Type Definitions */
/*@{ */
/** Pointer to char
 */
typedef const char* PointerToConstCharType;

/*@} */

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
//#define KRATOS_PRINT_TRACE
#ifdef KRATOS_PRINT_TRACE

#define KRATOS_TRACE(A,B) gTrace.Inform(A,B)

#else

#define KRATOS_TRACE(A,B)
#endif

#define KRATOS_WATCH(variable) \
  std::cout << #variable << " : " << variable << std::endl;

}  /* namespace Kratos.*/



#define KRATOS_SERIALIZE_SAVE_BASE_CLASS(Serializer, BaseType) \
	Serializer.save_base("BaseClass",*static_cast<const BaseType *>(this));

#define KRATOS_SERIALIZE_LOAD_BASE_CLASS(Serializer, BaseType) \
	Serializer.load_base("BaseClass",*static_cast<BaseType *>(this));


#endif /* KRATOS_DEFINE_H_INCLUDED  defined */
