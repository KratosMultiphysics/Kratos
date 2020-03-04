//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes


// External includes

// Project includes
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos
{

void FindNodalNeighboursProcess::Execute()
{
    mpNodeNeighboursCalculator->Execute();
    mpElemNeighboursCalculator->Execute();
}

/***********************************************************************************/
/***********************************************************************************/

void FindNodalNeighboursProcess::ClearNeighbours()
{
    mpNodeNeighboursCalculator->ClearNeighbours();
    mpElemNeighboursCalculator->ClearNeighbours();
}
  
}  // namespace Kratos.


