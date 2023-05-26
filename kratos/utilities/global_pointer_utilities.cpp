//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"
#include "utilities/global_pointer_utilities.h"

namespace Kratos
{
std::string GlobalPointerUtilities::Info() const
{
    std::stringstream buffer;
    buffer << "GlobalPointerUtilities" ;
    return buffer.str();
}

/***********************************************************************************/
/***********************************************************************************/

void GlobalPointerUtilities::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "GlobalPointerUtilities";
}

/***********************************************************************************/
/***********************************************************************************/

void GlobalPointerUtilities::PrintData(std::ostream& rOStream) const 
{

}

/***********************************************************************************/
/***********************************************************************************/

bool GlobalPointerUtilities::ObjectIsLocal(const Element& rElem, const int CurrentRank)
{
    return true; //if the iterator was found, then it is local!
}

/***********************************************************************************/
/***********************************************************************************/

bool GlobalPointerUtilities::ObjectIsLocal(const Condition& rCond, const int CurrentRank)
{
    return true; //if the iterator was found, then it is local!
}

/***********************************************************************************/
/***********************************************************************************/

bool GlobalPointerUtilities::ObjectIsLocal(const Node& rNode, const int CurrentRank)
{
    return rNode.FastGetSolutionStepValue(PARTITION_INDEX) == CurrentRank;
}

}  // namespace Kratos