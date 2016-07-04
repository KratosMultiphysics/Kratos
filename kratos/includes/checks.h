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

///@}

///@} addtogroup block

#endif // KRATOS_CHECKS_H_INCLUDED  defined
