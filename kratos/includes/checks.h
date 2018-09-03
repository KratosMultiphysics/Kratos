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
//

#include <cstring>
#include "includes/exception.h"

#if !defined(KRATOS_CHECKS_H_INCLUDED )
#define  KRATOS_CHECKS_H_INCLUDED

///@addtogroup KratosCore
///@{
///@name macros
///@{
#define KRATOS_CHECK(IsTrue) if(!(IsTrue)) KRATOS_ERROR << "Check failed because " << #IsTrue << " is not true"
#define KRATOS_CHECK_IS_FALSE(IsFalse) if(IsFalse) KRATOS_ERROR  << "Check failed because " << #IsFalse << " is not false"

#define KRATOS_CHECK_EQUAL(a,b) if(!(a == b)) KRATOS_ERROR << "Check failed because " << #a << " is not equal to " << #b
#define KRATOS_CHECK_NOT_EQUAL(a,b) if(a == b) KRATOS_ERROR << "Check failed because " << #a << " is equal to " << #b

#define KRATOS_CHECK_C_STRING_EQUAL(a,b) if((strcmp(a,b) != 0)) KRATOS_ERROR << "Check failed because \"" << a << "\" is not equal to \"" << b << "\""
#define KRATOS_CHECK_C_STRING_NOT_EQUAL(a,b) if((strcmp(a,b) == 0)) KRATOS_ERROR << "Check failed because \"" << a << "\" is equal to \"" << b << "\""

#define KRATOS_CHECK_LESS(a,b) if(!(a < b)) KRATOS_ERROR << "Check failed because " << #a << " is greater than or equal to " << #b
#define KRATOS_CHECK_LESS_EQUAL(a,b) if(a > b) KRATOS_ERROR << "Check failed because " << #a << " is greater than " << #b

#define KRATOS_CHECK_GREATER(a,b) if(!(a > b)) KRATOS_ERROR << "Check failed because " << #a << " is less than or equal to " << #b
#define KRATOS_CHECK_GREATER_EQUAL(a,b) if(a < b) KRATOS_ERROR  << "Check failed because " << #a << " is less than " << #b

#define KRATOS_CHECK_STRING_CONTAIN_SUB_STRING(TheString, SubString) if (TheString.find(SubString) == std::string::npos ) \
KRATOS_ERROR << "The string \"" << SubString << "\" was not found in the given string"

#define KRATOS_CHECK_NEAR(a,b, tolerance) if(std::abs(a - b) > tolerance) KRATOS_ERROR << "Check failed because " << #a << " = " << a << \
" is not near to " << #b << " = " << b << " within the tolerance " << tolerance
#define KRATOS_CHECK_DOUBLE_EQUAL(a,b) KRATOS_CHECK_NEAR(a,b,std::numeric_limits<double>::epsilon())

#define KRATOS_CHECK_EXCEPTION_IS_THROWN(TheStatement, TheErrorMessage)                 \
try {                                                                                   \
    TheStatement;                                                                       \
    KRATOS_ERROR << #TheStatement << " exited without throwing an error." << std::endl; \
} catch (Kratos::Exception& e) {                                                        \
    if ( std::string(e.what()).find( TheErrorMessage ) == std::string::npos )           \
        KRATOS_ERROR                                                                    \
            << "Test Failed: " << #TheStatement                                         \
            << " did not throw the expected error." << std::endl                        \
            << "Expected:" << std::endl << TheErrorMessage << std::endl                 \
            << "Got:" << std::endl << e.what() << std::endl;                            \
}

#define KRATOS_CHECK_VARIABLE_KEY(TheVariable)                               \
    KRATOS_ERROR_IF(TheVariable.Key() == 0)                                  \
        << TheVariable.Name() << " Key is 0." << std::endl                   \
        << "Check that Kratos variables have been correctly registered and " \
           "all required applications have been imported."                   \
        << std::endl;

#define KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TheVariable, TheNode)                          \
    KRATOS_ERROR_IF_NOT(TheNode.SolutionStepsDataHas(TheVariable))                         \
        << "Missing " << TheVariable.Name() << " variable in solution step data for node " \
        << TheNode.Id() << "." << std::endl;

#define KRATOS_CHECK_DOF_IN_NODE(TheVariable, TheNode)            \
    KRATOS_ERROR_IF_NOT(TheNode.HasDofFor(TheVariable))           \
        << "Missing Degree of Freedom for " << TheVariable.Name() \
        << " in node " << TheNode.Id() << "." << std::endl;

#ifdef KRATOS_DEBUG
#define KRATOS_DEBUG_CHECK(IsTrue) KRATOS_CHECK(IsTrue)
#define KRATOS_DEBUG_CHECK_IS_FALSE(IsFalse) KRATOS_CHECK_IS_FALSE(IsFalse)

#define KRATOS_DEBUG_CHECK_EQUAL(a,b) KRATOS_CHECK_EQUAL(a,b)
#define KRATOS_DEBUG_CHECK_NOT_EQUAL(a,b) KRATOS_CHECK_NOT_EQUAL(a,b)

#define KRATOS_DEBUG_CHECK_C_STRING_EQUAL(a,b) KRATOS_CHECK_C_STRING_EQUAL(a,b)
#define KRATOS_DEBUG_CHECK_C_STRING_NOT_EQUAL(a,b)  KRATOS_CHECK_C_STRING_NOT_EQUAL(a,b)

#define KRATOS_DEBUG_CHECK_LESS(a,b) KRATOS_CHECK_LESS(a,b)
#define KRATOS_DEBUG_CHECK_LESS_EQUAL(a,b) KRATOS_CHECK_LESS_EQUAL(a,b)

#define KRATOS_DEBUG_CHECK_GREATER(a,b) KRATOS_CHECK_GREATER(a,b)
#define KRATOS_DEBUG_CHECK_GREATER_EQUAL(a,b) KRATOS_CHECK_GREATER_EQUAL(a,b)

#define KRATOS_DEBUG_CHECK_STRING_CONTAIN_SUB_STRING(TheString, SubString) KRATOS_CHECK_STRING_CONTAIN_SUB_STRING(TheString, SubString)

#define KRATOS_DEBUG_CHECK_NEAR(a,b, tolerance) KRATOS_CHECK_NEAR(a,b, tolerance)
#define KRATOS_DEBUG_CHECK_DOUBLE_EQUAL(a,b) KRATOS_CHECK_DOUBLE_EQUAL(a,b)

#define KRATOS_DEBUG_CHECK_EXCEPTION_IS_THROWN(TheStatement, TheErrorMessage) KRATOS_CHECK_EXCEPTION_IS_THROWN(TheStatement, TheErrorMessage)

#define KRATOS_DEBUG_CHECK_VARIABLE_KEY(TheVariable) KRATOS_CHECK_VARIABLE_KEY(TheVariable)
#define KRATOS_DEBUG_CHECK_VARIABLE_IN_NODAL_DATA(TheVariable, TheNode) KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TheVariable, TheNode)
#define KRATOS_DEBUG_CHECK_DOF_IN_NODE(TheVariable, TheNode) KRATOS_CHECK_DOF_IN_NODE(TheVariable, TheNode)

#else
#define KRATOS_DEBUG_CHECK(IsTrue) if(false) KRATOS_CHECK(IsTrue)
#define KRATOS_DEBUG_CHECK_IS_FALSE(IsFalse) if(false) KRATOS_CHECK_IS_FALSE(IsFalse)

#define KRATOS_DEBUG_CHECK_EQUAL(a,b) if(false) KRATOS_CHECK_EQUAL(a,b)
#define KRATOS_DEBUG_CHECK_NOT_EQUAL(a,b) if(false) KRATOS_CHECK_NOT_EQUAL(a,b)

#define KRATOS_DEBUG_CHECK_C_STRING_EQUAL(a,b) if(false) KRATOS_CHECK_C_STRING_EQUAL(a,b)
#define KRATOS_DEBUG_CHECK_C_STRING_NOT_EQUAL(a,b) if(false)  KRATOS_CHECK_C_STRING_NOT_EQUAL(a,b)

#define KRATOS_DEBUG_CHECK_LESS(a,b) if(false) KRATOS_CHECK_LESS(a,b)
#define KRATOS_DEBUG_CHECK_LESS_EQUAL(a,b) if(false) KRATOS_CHECK_LESS_EQUAL(a,b)

#define KRATOS_DEBUG_CHECK_GREATER(a,b) if(false) KRATOS_CHECK_GREATER(a,b)
#define KRATOS_DEBUG_CHECK_GREATER_EQUAL(a,b) if(false) KRATOS_CHECK_GREATER_EQUAL(a,b)

#define KRATOS_DEBUG_CHECK_STRING_CONTAIN_SUB_STRING(TheString, SubString) if(false) KRATOS_CHECK_STRING_CONTAIN_SUB_STRING(TheString, SubString)

#define KRATOS_DEBUG_CHECK_NEAR(a,b, tolerance) if(false) KRATOS_CHECK_NEAR(a,b, tolerance)
#define KRATOS_DEBUG_CHECK_DOUBLE_EQUAL(a,b) if(false) KRATOS_CHECK_DOUBLE_EQUAL(a,b)

#define KRATOS_DEBUG_CHECK_EXCEPTION_IS_THROWN(TheStatement, TheErrorMessage) if(false) KRATOS_CHECK_EXCEPTION_IS_THROWN(TheStatement, TheErrorMessage)

#define KRATOS_DEBUG_CHECK_VARIABLE_KEY(TheVariable) if(false) KRATOS_CHECK_VARIABLE_KEY(TheVariable)
#define KRATOS_DEBUG_CHECK_VARIABLE_IN_NODAL_DATA(TheVariable, TheNode) if(false) KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TheVariable, TheNode)
#define KRATOS_DEBUG_CHECK_DOF_IN_NODE(TheVariable, TheNode) if(false) KRATOS_CHECK_DOF_IN_NODE(TheVariable, TheNode)
#endif
///@}

///@} addtogroup block

#endif // KRATOS_CHECKS_H_INCLUDED  defined
