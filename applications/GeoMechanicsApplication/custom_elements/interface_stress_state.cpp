// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//
#include "interface_stress_state.h"
#include "geo_mechanics_application_constants.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos
{

Matrix InterfaceStressState::CalculateBMatrix(const Matrix&, const Vector& rN, const Geometry<Node>& rGeometry) const
{
    KRATOS_ERROR_IF(rN.empty())
        << "Shape function values are empty, cannot compute the B matrix.\n";
    KRATOS_ERROR_IF_NOT(rN.size() == rGeometry.size() / 2)
        << "Each node pair must have exactly one shape function value, cannot compute the B "
           "matrix.\n";

    Matrix result = ZeroMatrix(GetVoigtSize(), rGeometry.WorkingSpaceDimension() * rGeometry.size());

    const auto number_of_u_dofs_per_side = result.size2() / 2;
    for (unsigned int i = 0; i < rGeometry.size() / 2; ++i) {
        result(0, i * rGeometry.WorkingSpaceDimension() + 1)                             = -rN[i];
        result(0, i * rGeometry.WorkingSpaceDimension() + 1 + number_of_u_dofs_per_side) = rN[i];

        result(1, i * rGeometry.WorkingSpaceDimension())                             = -rN[i];
        result(1, i * rGeometry.WorkingSpaceDimension() + number_of_u_dofs_per_side) = rN[i];
    }

    return result;
}

double InterfaceStressState::CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                             double DetJ,
                                                             const Geometry<Node>&) const
{
    return rIntegrationPoint.Weight() * DetJ;
}

Vector InterfaceStressState::CalculateGreenLagrangeStrain(const Matrix&) const
{
    KRATOS_ERROR << "For interfaces, it is not possible to calculate the Green Lagrange "
                    "strain based on a deformation gradient.\n";
}

std::unique_ptr<StressStatePolicy> InterfaceStressState::Clone() const
{
    return std::make_unique<InterfaceStressState>();
}

const Vector& InterfaceStressState::GetVoigtVector() const { return VoigtVectorInterface2D; }

SizeType InterfaceStressState::GetVoigtSize() const { return VOIGT_SIZE_2D_INTERFACE; }

SizeType InterfaceStressState::GetStressTensorSize() const { return 0; }

Vector InterfaceStressState::DefineInterfaceVoigtVector()
{
    Vector result{VOIGT_SIZE_2D_INTERFACE};
    result <<= 1.0, 0.0;

    return result;
}

const Vector InterfaceStressState::VoigtVectorInterface2D =
    InterfaceStressState::DefineInterfaceVoigtVector();
} // namespace Kratos
