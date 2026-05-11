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
 * @class PiecewiseLinearMomentCapacityConstitutiveLaw
 * @brief 1D moment-curvature law.
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) PiecewiseLinearMomentCapacityConstitutiveLaw : public ConstitutiveLaw
{
public:
    using BaseType = ConstitutiveLaw;

    KRATOS_CLASS_POINTER_DEFINITION(PiecewiseLinearMomentCapacityConstitutiveLaw);

    // Geometry constants
    static constexpr SizeType strain_size      = 3;
    static constexpr SizeType space_dimenstion = 3;

    PiecewiseLinearMomentCapacityConstitutiveLaw() = default;
    PiecewiseLinearMomentCapacityConstitutiveLaw(const PiecewiseLinearMomentCapacityConstitutiveLaw&) = default;
    PiecewiseLinearMomentCapacityConstitutiveLaw& operator=(const PiecewiseLinearMomentCapacityConstitutiveLaw&) = default;
    PiecewiseLinearMomentCapacityConstitutiveLaw(PiecewiseLinearMomentCapacityConstitutiveLaw&&) noexcept = default;
    PiecewiseLinearMomentCapacityConstitutiveLaw& operator=(PiecewiseLinearMomentCapacityConstitutiveLaw&&) noexcept = default;

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
