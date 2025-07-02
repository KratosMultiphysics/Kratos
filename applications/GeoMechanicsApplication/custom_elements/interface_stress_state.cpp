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
#include "includes/serializer.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos
{

Matrix Line2DInterfaceStressState::CalculateBMatrix(const Matrix&, const Vector& rN, const Geometry<Node>& rGeometry) const
{
    KRATOS_ERROR_IF(rN.empty())
        << "Shape function values are empty. Therefore, the B matrix can not be computed.\n";
    KRATOS_ERROR_IF_NOT(rN.size() == rGeometry.size() / 2)
        << "The number of shape functions should be equal to the number of node pairs. Therefore, "
           "the B matrix can not be computed.\n";

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

Vector Line2DInterfaceStressState::CalculateGreenLagrangeStrain(const Matrix&) const
{
    KRATOS_ERROR << "For interfaces, it is not possible to calculate the Green-Lagrange "
                    "strain based on a deformation gradient.\n";
}

std::unique_ptr<StressStatePolicy> Line2DInterfaceStressState::Clone() const
{
    return std::make_unique<Line2DInterfaceStressState>();
}

const Vector& Line2DInterfaceStressState::GetVoigtVector() const { return VoigtVectorInterface2D; }

SizeType Line2DInterfaceStressState::GetVoigtSize() const { return VOIGT_SIZE_2D_INTERFACE; }

SizeType Line2DInterfaceStressState::GetStressTensorSize() const
{
    KRATOS_ERROR << "For interfaces, the stress tensor size is not implemented.\n";
}

Vector Line2DInterfaceStressState::DefineInterfaceVoigtVector()
{
    Vector result{VOIGT_SIZE_2D_INTERFACE};
    result <<= 1.0, 0.0;

    return result;
}

const Vector Line2DInterfaceStressState::VoigtVectorInterface2D =
    Line2DInterfaceStressState::DefineInterfaceVoigtVector();

void Line2DInterfaceStressState::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void Line2DInterfaceStressState::load(Serializer&)
{
    // No data members to be loaded (yet)
}

} // namespace Kratos
