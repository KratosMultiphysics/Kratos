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
//                   Riccardo Rossi
//                    
//


#if !defined(KRATOS_GLOBAL_VARIABLES_H_INCLUDED )
#define  KRATOS_GLOBAL_VARIABLES_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "containers/variables_list.h"


namespace Kratos
{
namespace Globals
{

/*		class VariableData;
		class Element;
		class Condition;
	*/
///@name Kratos Globals
///@{
//
// This variable is NOT synchronized between different applications threads
extern KRATOS_API(KRATOS_CORE) VariablesList DefaultVariablesList;
/*
		extern KratosComponents<VariableData> VariableDataComponents;
		extern KratosComponents<Condition> ConditionComponents;
		extern KratosComponents<Element> ElementComponents;
  */
///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}

}

}  // namespace Kratos.

#endif // KRATOS_GLOBAL_VARIABLES_H_INCLUDED  defined
