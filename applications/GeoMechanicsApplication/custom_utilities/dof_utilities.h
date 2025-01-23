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

std::vector<std::size_t> KRATOS_API(GEO_MECHANICS_APPLICATION)
    ExtractEquationIdsFrom(const std::vector<Dof<double>*>& rDofs);

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
std::vector<Dof<double>*> ExtractDofsFromNodes(InputIt                 NodeRangeBegin,
                                               InputIt                 NodeRangeEnd,
                                               const Variable<double>& rDofVariable)
{
    auto result =
        std::vector<Dof<double>*>{static_cast<std::size_t>(std::distance(NodeRangeBegin, NodeRangeEnd))};
    ExtractDofsFromNodes(NodeRangeBegin, NodeRangeEnd, result.begin(), rDofVariable);
    return result;
}

template <typename NodeRange>
std::vector<Dof<double>*> ExtractDofsFromNodes(const NodeRange& rNodes, const Variable<double>& rDofVariable)
{
    return ExtractDofsFromNodes(std::begin(rNodes), std::end(rNodes), rDofVariable);
}

template <typename DisplacementNodeRange, typename WaterPressureNodeRange>
std::vector<Dof<double>*> ExtractUPwDofsFromNodes(const DisplacementNodeRange&  rDisplacementNodes,
                                                  const WaterPressureNodeRange& rWaterPressureNodes,
                                                  std::size_t                   ModelDimension)
{
    auto displacement_variables =
        Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y)};
    if (ModelDimension == 3) displacement_variables.push_back(std::cref(DISPLACEMENT_Z));

    auto result = std::vector<Dof<double>*>{
        displacement_variables.size() * std::size(rDisplacementNodes) + std::size(rWaterPressureNodes)};
    auto out_iter = result.begin();
    for (const auto& r_node : rDisplacementNodes) {
        out_iter = std::transform(
            std::begin(displacement_variables), std::end(displacement_variables), out_iter,
            [&r_node](const auto& r_variable) { return r_node.pGetDof(r_variable.get()); });
    }

    ExtractDofsFromNodes(std::begin(rWaterPressureNodes), std::end(rWaterPressureNodes), out_iter, WATER_PRESSURE);

    return result;
}

template <typename NodeRange>
std::vector<Dof<double>*> ExtractUPwDofsFromNodes(const NodeRange& rNodes, std::size_t ModelDimension)
{
    return ExtractUPwDofsFromNodes(rNodes, rNodes, ModelDimension);
}

Vector KRATOS_API(GEO_MECHANICS_APPLICATION)
    ExtractSolutionStepValues(const std::vector<Dof<double>*>& rDofs, int BufferIndex);
Vector KRATOS_API(GEO_MECHANICS_APPLICATION)
    ExtractFirstTimeDerivatives(const std::vector<Dof<double>*>& rDofs, int BufferIndex);
Vector KRATOS_API(GEO_MECHANICS_APPLICATION)
    ExtractSecondTimeDerivatives(const std::vector<Dof<double>*>& rDofs, int BufferIndex);

Vector KRATOS_API(GEO_MECHANICS_APPLICATION)
    ExtractSolutionStepValuesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int BufferIndex);
Vector KRATOS_API(GEO_MECHANICS_APPLICATION)
    ExtractFirstTimeDerivativesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int BufferIndex);
Vector KRATOS_API(GEO_MECHANICS_APPLICATION)
    ExtractSecondTimeDerivativesOfUPwDofs(const std::vector<Dof<double>*>& rDofs, int BufferIndex);

} // namespace Kratos::Geo::DofUtilities
