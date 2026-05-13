// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Gennady Markelov
//

#pragma once

// System includes

// External includes

#include "includes/constitutive_law.h"
#include "includes/table.h"

namespace Kratos
{

/**
 * @class PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw
 * @brief 1D moment-curvature law (plane-strain specific)
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw : public ConstitutiveLaw
{
public:
    using BaseType = ConstitutiveLaw;

    KRATOS_CLASS_POINTER_DEFINITION(PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw);

    // Geometry constants
    static constexpr SizeType strain_size     = 5;
    static constexpr SizeType space_dimension = 3;

    PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw() = default;
    PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw(
        const PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw&) = default;
    PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw& operator=(
        const PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw&) = default;
    PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw(
        PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw&&) noexcept = default;
    PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw& operator=(
        PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw&&) noexcept = default;

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    void GetLawFeatures(Features& rFeatures) override;

    double& CalculateValue(ConstitutiveLaw::Parameters& rParameters,
                           const Variable<double>&      rVariable,
                           double&                      rValue) override;
    using BaseType::CalculateValue;

    void CalculateMaterialResponsePK2(Parameters& rParameters) override;
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters) override;
    void FinalizeMaterialResponsePK2(Parameters& rParameters) override;

    void InitializeMaterial(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector&       rShapeFunctionsValues) override;

    [[nodiscard]] SizeType GetStrainSize() const override;

    bool RequiresFinalizeMaterialResponse() override;

    [[nodiscard]] int Check(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const ProcessInfo&  rCurrentProcessInfo) const override;

    std::string Info() const override;

private:
    Table<double, double> mStressStrainTable;
    double                mAccumulatedCurvature = 0.0;
    double                mUnReLoadCenter       = 0.0;
    double                mUnReLoadModulus      = 0.0;

    [[nodiscard]] double CalculateUnReLoadAmplitude() const;
    [[nodiscard]] bool   IsWithinUnReLoading(double Curvature) const;
    [[nodiscard]] std::pair<double, double> CalculateMomentAndTangentModulus(double curvature) const;
    [[nodiscard]] double BackboneTangentModulus(double curvature) const;
    [[nodiscard]] double BackboneMoment(double curvature) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos
