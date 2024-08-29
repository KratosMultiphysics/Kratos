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
    return {};
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

const Vector& InterfaceStressState::GetVoigtVector() const { return {}; }

SizeType InterfaceStressState::GetVoigtSize() const { return 0; }

SizeType InterfaceStressState::GetStressTensorSize() const { return 0; }

} // namespace Kratos
