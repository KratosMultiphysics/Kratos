//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes

// Application includes

// Include base h
#include "convection_diffusion_reaction_stabilization_adjoint_utilities.h"

namespace Kratos
{
namespace ConvectionDiffusionReactionStabilizationUtilities
{
template <>
double AdjointUtilities<2, 3>::CalculateElementLengthShapeDerivative(
    const double ElementLength,
    const double DetJDerivative)
{
    return 0.31830988618378353 * DetJDerivative / ElementLength;
}

template <>
double AdjointUtilities<3, 4>::CalculateElementLengthShapeDerivative(
    const double ElementLength,
    const double DetJDerivative)
{
    return 0.4714045207910277 * DetJDerivative / std::pow(ElementLength, 2);
}

// template instantiations
template class AdjointUtilities<2, 3>;
template class AdjointUtilities<3, 4>;

} // namespace ConvectionDiffusionReactionStabilizationUtilities
} // namespace Kratos
