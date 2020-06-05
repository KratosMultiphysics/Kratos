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


namespace Kratos
{
namespace Globals
{
	constexpr double Pi   = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651L;

    static constexpr int MaxAllowedThreads = 128; //we assume that no more than MaxAllowedThreads is used in SMP

/*		class VariableData;
		class Element;
		class Condition;
	*/
///@name Kratos Globals
///@{
//
// This variable is NOT synchronized between different applications threads
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
