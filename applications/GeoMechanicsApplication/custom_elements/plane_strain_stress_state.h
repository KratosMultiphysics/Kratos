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

#pragma once

#include "stress_state_policy.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) PlaneStrainStressState : public StressStatePolicy
{
public:
    [[nodiscard]] Matrix CalculateBMatrix(const Matrix&         rDN_DX,
                                          const Vector&         rN,
                                          const Geometry<Node>& rGeometry) const final;
    [[nodiscard]] double CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                         double DetJ,
                                                         const Geometry<Node>& rGeometry) const final;
    [[nodiscard]] Vector CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const final;
    [[nodiscard]] std::unique_ptr<StressStatePolicy> Clone() const final;
    [[nodiscard]] const Vector&                      GetVoigtVector() const final;
    [[nodiscard]] SizeType                           GetVoigtSize() const final;
    [[nodiscard]] SizeType                           GetStressTensorSize() const final;

private:
    [[nodiscard]] static Vector ConvertStrainTensorToVector(const Matrix& rStrainTensor);

    friend class Serializer;
    void save(Serializer&) const final;
    void load(Serializer&) final;
};

} // namespace Kratos
