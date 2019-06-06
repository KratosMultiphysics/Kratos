//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez

// Project includes
#include "custom_elements/incompressible_potential_flow_element.h"
#include "compressible_potential_flow_application_variables.h"

namespace Kratos
{
namespace PotentialFlowUtilities
{
template <int Dim, int NumNodes>
array_1d<double, NumNodes> GetPotentialOnNormalElement(const Element& rElement);

} // namespace PotentialFlow
} // namespace Kratos