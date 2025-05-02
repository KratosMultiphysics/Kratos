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

class KRATOS_API(GEO_MECHANICS_APPLICATION) PlaneStrainStressState final : public StressStatePolicy
{
public:
    [[nodiscard]] Matrix CalculateBMatrix(const Matrix&         rDN_DX,
                                          const Vector&         rN,
                                          const Geometry<Node>& rGeometry) const override;
    [[nodiscard]] Vector CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const override;
    [[nodiscard]] std::unique_ptr<StressStatePolicy> Clone() const override;
    [[nodiscard]] const Vector&                      GetVoigtVector() const override;
    [[nodiscard]] SizeType                           GetVoigtSize() const override;
    [[nodiscard]] SizeType                           GetStressTensorSize() const override;

private:
    [[nodiscard]] static Vector ConvertStrainTensorToVector(const Matrix& rStrainTensor);

    friend class Serializer;
    void save(Serializer&) const override;
    void load(Serializer&) override;
};

} // namespace Kratos
