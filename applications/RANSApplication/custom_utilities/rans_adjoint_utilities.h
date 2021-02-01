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

#if !defined(KRATOS_RANS_APPLICATION_ADJOINT_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_ADJOINT_UTILITIES_H_INCLUDED

// System includes
#include <cmath>
#include <tuple>

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/model_part.h"
#include "utilities/geometrical_sensitivity_utility.h"

namespace Kratos
{
///@name Kratos Globals
///@{

class AdjointUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using ConditionType = ModelPart::ConditionType;

    ///@}
    ///@name Static Operations
    ///@{

    static double CalculateWallHeightDerivative(
        const ConditionType& rCondition,
        const array_1d<double, 3>& rNormal,
        const array_1d<double, 3>& rNormalDerivative);

    static void CalculateYPlusAndUtauDerivative(
        double& rYPlusDerivative,
        double& rUTauDerivative,
        const double rYPlus,
        const double rUTau,
        const double WallVelocity,
        const double WallVelocityDerivative,
        const double WallHeight,
        const double WallHeightDerivative,
        const double KinematicViscosity,
        const double Kappa,
        const double Beta,
        const double YPlusLimit);

    ///@}
};

} // namespace Kratos

#endif // KRATOS_RANS_APPLICATION_ADJOINT_UTILITIES_H_INCLUDED