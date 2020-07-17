//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes
#include <cmath>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"

// Application includes

// Include base h
#include "convection_diffusion_reaction_stabilization_adjoint_utilities.h"

namespace Kratos
{
///@name  Functions
///@{

namespace ConvectionDiffusionReactionStabilizationUtilities
{
template <>
double AdjointUtilities<2, 3>::CalculateElementLengthShapeDerivative(const double ElementLength,
                                                                     const double DetJDerivative)
{
    return 0.31830988618378353 * DetJDerivative / ElementLength;
}

template <>
double AdjointUtilities<3, 4>::CalculateElementLengthShapeDerivative(const double ElementLength,
                                                                     const double DetJDerivative)
{
    return 0.4714045207910277 * DetJDerivative / std::pow(ElementLength, 2);
}

} // namespace ConvectionDiffusionReactionStabilizationUtilities
} // namespace Kratos