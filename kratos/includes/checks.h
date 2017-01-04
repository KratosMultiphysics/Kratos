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
#define KRATOS_CHECK(IsTrue) if(!(IsTrue)) KRATOS_ERROR << "Check failed becuase " << #IsTrue << " is not true"
#define KRATOS_CHECK_IS_FALSE(IsFalse) if(IsFalse) KRATOS_ERROR  << "Check failed becuase " << #IsFalse << " is not false"

#define KRATOS_CHECK_EQUAL(a,b) if(!(a == b)) KRATOS_ERROR << "Check failed becuase " << #a << " is not equal to " << #b
#define KRATOS_CHECK_NOT_EQUAL(a,b) if(a == b) KRATOS_ERROR << "Check failed becuase " << #a << " is equal to " << #b

#define KRATOS_CHECK_C_STRING_EQUAL(a,b) if((strcmp(a,b) != 0)) KRATOS_ERROR << "Check failed becuase \"" << a << "\" is not equal to \"" << b << "\""
#define KRATOS_CHECK_C_STRING_NOT_EQUAL(a,b) if((strcmp(a,b) == 0)) KRATOS_ERROR << "Check failed becuase \"" << a << "\" is equal to \"" << b << "\""

#define KRATOS_CHECK_LESS(a,b) if(!(a < b)) KRATOS_ERROR << "Check failed becuase " << #a << " is greater than or equal to " << #b
#define KRATOS_CHECK_LESS_EQUAL(a,b) if(a > b) KRATOS_ERROR << "Check failed becuase " << #a << " is greater than " << #b

#define KRATOS_CHECK_GREATER(a,b) if(!(a > b)) KRATOS_ERROR << "Check failed becuase " << #a << " is less than or equal to " << #b
#define KRATOS_CHECK_GREATER_EQUAL(a,b) if(a < b) KRATOS_ERROR  << "Check failed becuase " << #a << " is less than " << #b

#define KRATOS_CHECK_STRING_CONTAIN_SUB_STRING(TheString, SubString) if (TheString.find(SubString) == std::string::npos ) \
KRATOS_ERROR << "The string \"" << SubString << "\" was not found in the given string"

#define KRATOS_CHECK_NEAR(a,b, tolerance) if(std::abs(a - b) > tolerance) KRATOS_ERROR << "Check failed becuase " << #a << " = " << a << \
" is not near to " << #b << " = " << b << " within the tolerance " << tolerance
#define KRATOS_CHECK_DOUBLE_EQUAL(a,b) KRATOS_CHECK_NEAR(a,b,std::numeric_limits<double>::epsilon())

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
#endif
///@}

///@} addtogroup block

#endif // KRATOS_CHECKS_H_INCLUDED  defined
