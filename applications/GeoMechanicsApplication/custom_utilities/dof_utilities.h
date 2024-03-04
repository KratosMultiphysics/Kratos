// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#pragma once

#include <algorithm>
#include <vector>

#include "geometries/geometry.h"
#include "includes/dof.h"
#include "includes/node.h"
#include "includes/variables.h"

namespace Kratos::Geo::DofUtilities
{

std::vector<std::size_t> ExtractEquationIdsFrom(const std::vector<Dof<double>*>& rDofs);

template <typename InputIt>
std::vector<Dof<double>*> ExtractDofsFromNodes(InputIt                 NodePtrRangeBegin,
                                               InputIt                 NodePtrRangeEnd,
                                               const Variable<double>& rDofVariable)
{
    auto result = std::vector<Dof<double>*>{};
    std::transform(NodePtrRangeBegin, NodePtrRangeEnd, std::back_inserter(result),
                   [&rDofVariable](const auto& r_node) { return r_node.pGetDof(rDofVariable); });
    return result;
}

template <typename NodePtrRange>
std::vector<Dof<double>*> ExtractDofsFromNodes(const NodePtrRange& rNodePtrs, const Variable<double>& rDofVariable)
{
    return ExtractDofsFromNodes(std::begin(rNodePtrs), std::end(rNodePtrs), rDofVariable);
}

template <typename NodeRange>
std::vector<Dof<double>*> ExtractUPwDofsFromNodes(const NodeRange& rNodes)
{
    auto result = std::vector<Dof<double>*>{};
    for (auto& r_node : rNodes) {
        result.push_back(r_node.pGetDof(DISPLACEMENT_X));
        result.push_back(r_node.pGetDof(DISPLACEMENT_Y));
    }

    for (auto& r_node : rNodes) {
        result.push_back(r_node.pGetDof(WATER_PRESSURE));
    }

    return result;
}

Vector ExtractSolutionStepValues(const std::vector<Dof<double>*>& rDofs, int BufferIndex);
Vector ExtractFirstTimeDerivatives(const std::vector<Dof<double>*>& rDofs, int BufferIndex);
Vector ExtractSecondTimeDerivatives(const std::vector<Dof<double>*>& rDofs, int BufferIndex);

Vector ExtractSolutionStepValuesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int BufferIndex);
Vector ExtractFirstTimeDerivativesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int BufferIndex);
Vector ExtractSecondTimeDerivativesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int BufferIndex);

} // namespace Kratos::Geo::DofUtilities
