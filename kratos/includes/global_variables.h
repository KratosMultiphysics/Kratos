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
	constexpr double Pi   = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651L;
	constexpr double Pi_2 = 1.57079632679489661923132169163975144209858469968755291048747229615390820314310449931401741267105853399107404326L;
	constexpr double Pi_3 = 1.04719755119659774615421446109316762806572313312503527365831486410260546876206966620934494178070568932738269550L;

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
