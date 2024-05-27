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

// System includes
#include <cmath>
#include <limits>
#include <cstring>

// External includes

// Project includes
#include "includes/exception.h"

#pragma once

///@addtogroup KratosCore
///@{
///@name macros
///@{
#define KRATOS_EXPECT_TRUE(IsTrue) if(!(IsTrue)) KRATOS_ERROR << "Check failed because " << #IsTrue << " is not true" << std::endl;
#define KRATOS_EXPECT_FALSE(IsFalse) if(IsFalse) KRATOS_ERROR  << "Check failed because " << #IsFalse << " is not false" << std::endl;

#define KRATOS_EXPECT_EQ(a,b) if(!((a) == (b))) KRATOS_ERROR << "Check failed because " << #a << " is not equal to " << #b
#define KRATOS_EXPECT_NE(a,b) if((a) == (b)) KRATOS_ERROR << "Check failed because " << #a << " is equal to " << #b

// #define KRATOS_EXPECT_STREQ(a,b) if(a.compare(b) != 0) KRATOS_ERROR << "Check failed because \"" << a << "\" is not equal to \"" << b << "\"" << std::endl;
// #define KRATOS_EXPECT_STRNE(a,b) if(a.compare(b) == 0) KRATOS_ERROR << "Check failed because \"" << a << "\" is equal to \"" << b << "\"" << std::endl;

#define KRATOS_EXPECT_STREQ(a,b) if((strcmp(a,b) != 0)) KRATOS_ERROR << "Check failed because \"" << a << "\" is not equal to \"" << b << "\"" << std::endl;
#define KRATOS_EXPECT_STRNE(a,b) if((strcmp(a,b) == 0)) KRATOS_ERROR << "Check failed because \"" << a << "\" is equal to \"" << b << "\"" << std::endl;

#define KRATOS_EXPECT_LT(a,b) if(!(a < b)) KRATOS_ERROR << "Check failed because " << #a << " is greater than or equal to " << #b << std::endl;
#define KRATOS_EXPECT_LE(a,b) if(!(a <= b)) KRATOS_ERROR << "Check failed because " << #a << " is greater than " << #b << std::endl;

#define KRATOS_EXPECT_GT(a,b) if(!(a > b)) KRATOS_ERROR << "Check failed because " << #a << " is less than or equal to " << #b
#define KRATOS_EXPECT_GE(a,b) if(!(a >= b)) KRATOS_ERROR  << "Check failed because " << #a << " is less than " << #b

#define KRATOS_EXPECT_HAS_SUBSTRING(TheString, SubString) if (TheString.find(SubString) == std::string::npos ) \
KRATOS_ERROR << "The string \"" << SubString << "\" was not found in the given string" << std::endl;

#define KRATOS_EXPECT_NEAR(a,b, tolerance) if(!(std::abs(a - b) <= tolerance)) KRATOS_ERROR << "Check failed because " << #a << " = " << a << \
" is not near to " << #b << " = " << b << " within the tolerance " << tolerance
#define KRATOS_EXPECT_RELATIVE_NEAR(a,b, tolerance) if(!(std::abs(b) <= std::numeric_limits<double>::epsilon())) { KRATOS_ERROR_IF(!(std::abs((a - b)/b) <= tolerance)) << "Check failed because " << #a << " = " << a << \
" is not near to " << #b << " = " << b << " within the relative tolerance " << tolerance << std::endl; } else {KRATOS_EXPECT_NEAR(a,b,tolerance);}
#define KRATOS_EXPECT_DOUBLE_EQ(a,b) KRATOS_EXPECT_NEAR(a,b,std::numeric_limits<double>::epsilon())

#define KRATOS_EXPECT_VECTOR_NEAR(a, b, tolerance) {                        \
KRATOS_ERROR_IF_NOT(a.size() == b.size())                                  \
<< "Check failed because vector arguments do not have the same size:"      \
<< std::endl                                                               \
<< "First argument has size " << a.size() << ", "                          \
<< "second argument has size " << b.size() << "." << std::endl;            \
for (std::size_t _i = 0; _i < a.size(); _i++) {                              \
   KRATOS_ERROR_IF( !(std::abs(a[_i] - b[_i]) <= tolerance) )                \
   << "Check failed because vector " << #a << " with values" << std::endl  \
   << a << std::endl                                                       \
   << "Is not near vector " << #b << " with values" << std::endl           \
   << b << std::endl                                                       \
   << "Mismatch found in component " << _i << ":" << std::endl              \
   << a[_i] << " not near " << b[_i]                                         \
   << " within tolerance " << tolerance << "." << std::endl;               \
}                                                                          \
}

#define KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(a, b, tolerance) {                   \
KRATOS_ERROR_IF_NOT(a.size() == b.size())                                      \
<< "Check failed because vector arguments do not have the same size:"          \
<< std::endl                                                                   \
<< "First argument has size " << a.size() << ", "                              \
<< "second argument has size " << b.size() << "." << std::endl;                \
for (std::size_t _i = 0; _i < a.size(); _i++) {                                   \
    if (std::abs(b[_i]) > std::numeric_limits<double>::epsilon()) {             \
       KRATOS_ERROR_IF( !(std::abs((a[_i] - b[_i])/b[_i]) <= tolerance) )         \
       << "Check failed because vector " << #a << " with values" << std::endl  \
       << a << std::endl                                                       \
       << "Is not near vector " << #b << " with values" << std::endl           \
       << b << std::endl                                                       \
       << "Mismatch found in component " << _i << ":" << std::endl              \
       << a[_i] << " not near " << b[_i]                                         \
       << " within relative tolerance " << tolerance << "." << std::endl;      \
    } else {                                                                   \
       KRATOS_ERROR_IF( !(std::abs(a[_i] - b[_i]) <= tolerance) )                \
       << "Check failed because vector " << #a << " with values" << std::endl  \
       << a << std::endl                                                       \
       << "Is not near vector " << #b << " with values" << std::endl           \
       << b << std::endl                                                       \
       << "Mismatch found in component " << _i << ":" << std::endl              \
       << a[_i] << " not near " << b[_i]                                         \
       << " within tolerance " << tolerance << "." << std::endl;               \
    }                                                                          \
}                                                                              \
}

#define KRATOS_EXPECT_VECTOR_EQ(a, b) KRATOS_EXPECT_VECTOR_NEAR(a,b,std::numeric_limits<double>::epsilon())

#define KRATOS_EXPECT_MATRIX_NEAR(a, b, tolerance) {                              \
KRATOS_ERROR_IF_NOT((a.size1() == b.size1()) && (a.size2() == b.size2()))        \
<< "Check failed because matrix arguments do not have the same dimensions:"      \
<< std::endl                                                                     \
<< "First argument has dimensions (" << a.size1() << "," << a.size2() << "), "   \
<< "second argument has dimensions (" << b.size1() << "," << b.size2() << ")."   \
<< std::endl;                                                                    \
for (std::size_t _i = 0; _i < a.size1(); _i++) {                                    \
    for (std::size_t _j = 0; _j < a.size2(); _j++) {                                \
       KRATOS_ERROR_IF( !(std::abs(a(_i,_j) - b(_i,_j)) <= tolerance) )              \
       << "Check failed because matrix " << #a << " with values" << std::endl    \
       << a << std::endl                                                         \
       << "Is not near matrix " << #b << " with values" << std::endl             \
       << b << std::endl                                                         \
       << "Mismatch found in component (" << _i << "," << _j << "): " << std::endl \
       << a(_i,_j) << " not near " << b(_i,_j)                                       \
       << " within tolerance " << tolerance << "." << std::endl;                 \
    }                                                                            \
}                                                                                \
}

#define KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(a, b, tolerance) {                         \
KRATOS_ERROR_IF_NOT((a.size1() == b.size1()) && (a.size2() == b.size2()))            \
<< "Check failed because matrix arguments do not have the same dimensions:"          \
<< std::endl                                                                         \
<< "First argument has dimensions (" << a.size1() << "," << a.size2() << "), "       \
<< "second argument has dimensions (" << b.size1() << "," << b.size2() << ")."       \
<< std::endl;                                                                        \
for (std::size_t _i = 0; _i < a.size1(); _i++) {                                        \
    for (std::size_t _j = 0; _j < a.size2(); _j++) {                                    \
        if (std::abs(b(_i,_j)) > std::numeric_limits<double>::epsilon()) {             \
           KRATOS_ERROR_IF( !(std::abs((a(_i,_j) - b(_i,_j))/b(_i,_j)) <= tolerance) )     \
           << "Check failed because matrix " << #a << " with values" << std::endl    \
           << a << std::endl                                                         \
           << "Is not near matrix " << #b << " with values" << std::endl             \
           << b << std::endl                                                         \
           << "Mismatch found in component (" << _i << "," << _j << "): " << std::endl \
           << a(_i,_j) << " not near " << b(_i,_j)                                       \
           << " within relative tolerance " << tolerance << "." << std::endl;        \
        } else {                                                                     \
           KRATOS_ERROR_IF( !(std::abs(a(_i,_j) - b(_i,_j)) <= tolerance) )              \
           << "Check failed because matrix " << #a << " with values" << std::endl    \
           << a << std::endl                                                         \
           << "Is not near matrix " << #b << " with values" << std::endl             \
           << b << std::endl                                                         \
           << "Mismatch found in component (" << _i << "," << _j << "): " << std::endl \
           << a(_i,_j) << " not near " << b(_i,_j)                                       \
           << " within tolerance " << tolerance << "." << std::endl;                 \
    }                                                                                \
}                                                                                    \
}                                                                                    \
}

#define KRATOS_EXPECT_MATRIX_EQ(a, b) KRATOS_EXPECT_MATRIX_NEAR(a,b,std::numeric_limits<double>::epsilon())

#define KRATOS_EXPECT_EXCEPTION_IS_THROWN(TheStatement, TheErrorMessage)                 \
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

#define KRATOS_EXPECT_VARIABLE_IN_NODAL_DATA(TheVariable, TheNode)                          \
    KRATOS_ERROR_IF_NOT(TheNode.SolutionStepsDataHas(TheVariable))                         \
        << "Missing " << TheVariable.Name() << " variable in solution step data for node " \
        << TheNode.Id() << "." << std::endl;

#define KRATOS_EXPECT_DOF_IN_NODE(TheVariable, TheNode)            \
    KRATOS_ERROR_IF_NOT(TheNode.HasDofFor(TheVariable))           \
        << "Missing Degree of Freedom for " << TheVariable.Name() \
        << " in node " << TheNode.Id() << "." << std::endl;

#ifdef KRATOS_DEBUG
#define KRATOS_DEBUG_EXCEPT_EXCEPTION_IS_THROWN(TheStatement, TheErrorMessage) KRATOS_EXPECT_EXCEPTION_IS_THROWN(TheStatement, TheErrorMessage)
#else
#define KRATOS_DEBUG_EXCEPT_EXCEPTION_IS_THROWN(TheStatement, TheErrorMessage) if(false) KRATOS_EXPECT_EXCEPTION_IS_THROWN(TheStatement, TheErrorMessage)
#endif

///@}
///@} addtogroup block
