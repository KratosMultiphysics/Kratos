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
//                   Marjan Fathian
//

#pragma once

#include "stress_state_policy.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) AxisymmetricStressState : public StressStatePolicy
{
public:
    [[nodiscard]] double CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                         double DetJ,
                                                         const Geometry<Node>& rGeometry) const final;
    [[nodiscard]] Matrix CalculateBMatrix(const Matrix&         rDN_DX,
                                          const Vector&         rN,
                                          const Geometry<Node>& rGeometry) const final;
    [[nodiscard]] Vector CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const final;
    [[nodiscard]] std::unique_ptr<StressStatePolicy> Clone() const final;
    [[nodiscard]] const Vector&                      GetVoigtVector() const final;
    [[nodiscard]] SizeType                           GetVoigtSize() const final;
    [[nodiscard]] SizeType                           GetStressTensorSize() const final;

private:
    friend class Serializer;
    void save(Serializer&) const final;
    void load(Serializer&) final;
};

} // namespace Kratos
