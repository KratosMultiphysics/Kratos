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

#include "geo_aliases.h"
#include "geometries/geometry.h"
#include "includes/dof.h"
#include "includes/node.h"
#include "includes/variables.h"

namespace Kratos::Geo::DofUtilities
{

std::vector<std::size_t> ExtractEquationIdsFrom(const std::vector<Dof<double>*>& rDofs);

template <typename InputIt, typename OutputIt>
OutputIt ExtractDofsFromNodes(InputIt                 NodeRangeBegin,
                              InputIt                 NodeRangeEnd,
                              OutputIt                DofPtrRangeBegin,
                              const Variable<double>& rDofVariable)
{
    return std::transform(NodeRangeBegin, NodeRangeEnd, DofPtrRangeBegin, [&rDofVariable](const auto& r_node) {
        return r_node.pGetDof(rDofVariable);
    });
}

template <typename InputIt>
std::vector<Dof<double>*> ExtractDofsFromNodes(InputIt                 NodePtrRangeBegin,
                                               InputIt                 NodePtrRangeEnd,
                                               const Variable<double>& rDofVariable)
{
    auto result = std::vector<Dof<double>*>{};
    ExtractDofsFromNodes(NodePtrRangeBegin, NodePtrRangeEnd, std::back_inserter(result), rDofVariable);
    return result;
}

template <typename NodePtrRange>
std::vector<Dof<double>*> ExtractDofsFromNodes(const NodePtrRange& rNodePtrs, const Variable<double>& rDofVariable)
{
    return ExtractDofsFromNodes(std::begin(rNodePtrs), std::end(rNodePtrs), rDofVariable);
}

template <typename NodeRange1, typename NodeRange2>
std::vector<Dof<double>*> ExtractUPwDofsFromNodes(const NodeRange2& rSecondOrderNodes,
                                                  const NodeRange1& rFirstOrderNodes,
                                                  std::size_t       ModelDimension)
{
    auto result = std::vector<Dof<double>*>{};
    auto displacement_variables =
        Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y)};
    if (ModelDimension == 3) displacement_variables.push_back(std::cref(DISPLACEMENT_Z));

    for (const auto& r_node : rSecondOrderNodes) {
        for (const auto& r_variable : displacement_variables) {
            result.push_back(r_node.pGetDof(r_variable.get()));
        }
    }

    ExtractDofsFromNodes(std::begin(rFirstOrderNodes), std::end(rFirstOrderNodes),
                         std::back_inserter(result), WATER_PRESSURE);

    return result;
}

template <typename NodeRange>
std::vector<Dof<double>*> ExtractUPwDofsFromNodes(const NodeRange& rNodes, std::size_t ModelDimension)
{
    return ExtractUPwDofsFromNodes(rNodes, rNodes, ModelDimension);
}

Vector ExtractSolutionStepValues(const std::vector<Dof<double>*>& rDofs, int BufferIndex);
Vector ExtractFirstTimeDerivatives(const std::vector<Dof<double>*>& rDofs, int BufferIndex);
Vector ExtractSecondTimeDerivatives(const std::vector<Dof<double>*>& rDofs, int BufferIndex);

Vector ExtractSolutionStepValuesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int BufferIndex);
Vector ExtractFirstTimeDerivativesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int BufferIndex);
Vector ExtractSecondTimeDerivativesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int BufferIndex);

} // namespace Kratos::Geo::DofUtilities
