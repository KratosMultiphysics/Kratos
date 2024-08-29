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

namespace Kratos
{

Matrix InterfaceStressState::CalculateBMatrix(const Matrix& rDN_DX, const Vector& rN, const Geometry<Node>& rGeometry) const
{
    if (rN.empty()) return {};

    Matrix result = ZeroMatrix(GetVoigtSize(), rGeometry.WorkingSpaceDimension() * rGeometry.size());

    for (unsigned int i = 0; i < rGeometry.size() / 2; ++i) {
        result(0, i * rGeometry.WorkingSpaceDimension())     = 0;
        result(0, i * rGeometry.WorkingSpaceDimension() + 1) = -rN[i];

        result(1, i * rGeometry.WorkingSpaceDimension())     = -rN[i];
        result(1, i * rGeometry.WorkingSpaceDimension() + 1) = 0;
    }

    for (unsigned int i = rGeometry.size() / 2; i < rGeometry.size(); ++i) {
        result(0, i * rGeometry.WorkingSpaceDimension())     = 0;
        result(0, i * rGeometry.WorkingSpaceDimension() + 1) = rN[i - rGeometry.size() / 2];

        result(1, i * rGeometry.WorkingSpaceDimension())     = rN[i - rGeometry.size() / 2];
        result(1, i * rGeometry.WorkingSpaceDimension() + 1) = 0;
    }

    return result;
}

double InterfaceStressState::CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                             double                DetJ,
                                                             const Geometry<Node>& rGeometry) const
{
    return 0.0;
}

Vector InterfaceStressState::CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const
{
    return {};
}

std::unique_ptr<StressStatePolicy> InterfaceStressState::Clone() const
{
    return std::make_unique<InterfaceStressState>();
}

const Vector& InterfaceStressState::GetVoigtVector() const { return VoigtVectorInterface2D; }

SizeType InterfaceStressState::GetVoigtSize() const { return GetVoigtSizeInterface2D(); }

SizeType InterfaceStressState::GetStressTensorSize() const { return 0; }

} // namespace Kratos
