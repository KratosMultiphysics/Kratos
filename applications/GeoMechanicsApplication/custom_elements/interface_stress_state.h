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
//                   Anne van de Graaf
//

#pragma once

#include "stress_state_policy.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) Line2DInterfaceStressState final : public StressStatePolicy
{
public:
    [[nodiscard]] Matrix CalculateBMatrix(const Matrix&, const Vector& rN, const Geometry<Node>& rGeometry) const override;
    [[noreturn]] Vector CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const override;
    [[nodiscard]] std::unique_ptr<StressStatePolicy> Clone() const override;
    [[nodiscard]] const Vector&                      GetVoigtVector() const override;
    [[nodiscard]] SizeType                           GetVoigtSize() const override;
    [[noreturn]] SizeType                            GetStressTensorSize() const override;

private:
    static const Vector VoigtVectorInterface2D;
    static Vector       DefineInterfaceVoigtVector();

    friend class Serializer;
    void save(Serializer&) const override;
    void load(Serializer&) override;
};

class KRATOS_API(GEO_MECHANICS_APPLICATION) SurfaceInterfaceStressState final : public StressStatePolicy
{
public:
    [[nodiscard]] Matrix CalculateBMatrix(const Matrix&, const Vector& rN, const Geometry<Node>& rGeometry) const override;
    [[noreturn]] Vector CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const override;
    [[nodiscard]] std::unique_ptr<StressStatePolicy> Clone() const override;
    [[nodiscard]] const Vector&                      GetVoigtVector() const override;
    [[nodiscard]] SizeType                           GetVoigtSize() const override;
    [[noreturn]] SizeType                            GetStressTensorSize() const override;

private:
    static const Vector VoigtVectorInterface3D;
    static Vector       DefineInterfaceVoigtVector();

    friend class Serializer;
    void save(Serializer&) const override;
    void load(Serializer&) override;
};

} // namespace Kratos
