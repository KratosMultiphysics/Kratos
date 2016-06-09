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

#include "includes/kratos_exception.h"

#if !defined(KRATOS_CHECKS_H_INCLUDED )
#define  KRATOS_CHECKS_H_INCLUDED

///@addtogroup KratosCore
///@{
///@name macros
///@{ 
#define KRATOS_CHECK(IsTrue) if(!(IsTrue)) KRATOS_ERROR 
#define KRATOS_CHECK_IS_FALSE(IsFalse) if(IsFalse) KRATOS_ERROR 

#define KRATOS_CHECK_EQUAL(a,b) if(!(a == b)) KRATOS_ERROR 
#define KRATOS_CHECK_NOT_EQUAL(a,b) if(a == b) KRATOS_ERROR 

#define KRATOS_CHECK_LESS(a,b) if(!(a < b)) KRATOS_ERROR 
#define KRATOS_CHECK_LESS_EQUAL(a,b) if(a > b) KRATOS_ERROR 

#define KRATOS_CHECK_GREATER(a,b) if(!(a > b)) KRATOS_ERROR 
#define KRATOS_CHECK_GREATER_EQUAL(a,b) if(a < b) KRATOS_ERROR 


///@}

///@} addtogroup block

#endif // KRATOS_CHECKS_H_INCLUDED  defined 


