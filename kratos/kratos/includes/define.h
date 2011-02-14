/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 

/* *********************************************************   
 *          
 *   Last modified by:    $Author: rrossi $
 *   Date:                $Date: 2007-03-06 10:30:33 $
 *   Revision:            $Revision: 1.2 $
 *
 * ***********************************************************/


#if !defined(KRATOS_DEFINE )
#define  KRATOS_DEFINE


/* System includes */
#include <stdexcept>
#include <sstream>


/* External includes */
#include "boost/smart_ptr.hpp"
#include "boost/current_function.hpp"


/* Project includes */

#include "includes/kratos_config.h"


#define KRATOS_CLASS_POINTER_DEFINITION(a) typedef boost::shared_ptr<a > Pointer; \
typedef boost::shared_ptr<a > SharedPointer; \
typedef boost::weak_ptr<a > WeakPointer


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

#define KRATOS_CATCH_AND_THROW(ExceptionType, MoreInfo, Block) \
catch(ExceptionType& e)                                        \
{                                                              \
Block                                                          \
std::stringstream buffer;                                      \
buffer << e.what() << std::endl;                               \
buffer << "\nwhile executing : " << BOOST_CURRENT_FUNCTION << " [ " << __FILE__ << " , Line " << __LINE__ << " ] " << MoreInfo ; \
throw ExceptionType(buffer.str());                             \
}

#define KRATOS_ERROR(ExceptionType, ErrorMessage, MoreInfo)    \
{                                                              \
std::stringstream kratos_error_buffer_12345;                                      \
kratos_error_buffer_12345 << "in " << BOOST_CURRENT_FUNCTION << " [ " << __FILE__ << " , Line " << __LINE__ << " ]" << std::endl; \
kratos_error_buffer_12345 << "\nwith subject    :  " << ErrorMessage << " " << MoreInfo; \
throw ExceptionType(kratos_error_buffer_12345.str());                             \
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
catch(std::exception& e) { Block KRATOS_ERROR(std::runtime_error, e.what(), MoreInfo) } \
catch(...) { Block KRATOS_ERROR(std::runtime_error, "Unknown error", MoreInfo) } 

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

#ifdef KRATOS_DEFINE_VARIABLE
#undef KRATOS_DEFINE_VARIABLE
#endif
#define KRATOS_DEFINE_VARIABLE(type, name) \
    extern /*const*/ Variable<type > name;

#ifdef KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS
#undef KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(name) \
    extern /*const*/ Kratos::Variable<Kratos::array_1d<double, 3> > name; \
    extern /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_X;\
    extern /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_Y;\
    extern /*const*/ Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > name##_Z;

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
    AddComponent(name.Name(), name); \
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
// components
//
//-----------------------------------------------------------------

#ifdef KRATOS_REGISTER_ELEMENT
#undef KRATOS_REGISTER_ELEMENT
#endif
#define KRATOS_REGISTER_ELEMENT(name, reference) \
    KratosComponents<Element >::Add(name, reference);

#ifdef KRATOS_REGISTER_CONDITION
#undef KRATOS_REGISTER_CONDITION
#endif
#define KRATOS_REGISTER_CONDITION(name, reference) \
    KratosComponents<Condition >::Add(name, reference);





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

#define KRATOS_TRIANGULAR_MEMBRANE_ELEMENT_INCLUDED
#define KRATOS_QUADRILATERAL_DIFFUSION_CONVECTION_ELEMENT_INCLUDED
#define KRATOS_TETRAHEDRAL_HEAT_CONDUCTIVITY_ELEMENT_INCLUDED

#define KRATOS_WATCH(variable) \
  std::cout << #variable << " : " << variable << std::endl;

  /* Major version
   */
#define KRATOS_MAJOR_VERSION  0

  /* Minor version 
   */
#define KRATOS_MINOR_VERSION  2

  /* Corrected version 
   */
#define KRATOS_CORRECTION     0
  

}  /* namespace Kratos.*/

#define KRATOS_SERIALIZE_SAVE_BASE_CLASS(Serializer, BaseType)

#define KRATOS_SERIALIZE_LOAD_BASE_CLASS(Serializer, BaseType)

  
#endif /* KRATOS_DEFINE  defined */

