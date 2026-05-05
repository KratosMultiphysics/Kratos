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
 * @brief 1D moment-curvature law with four regimes and axial-load dependent capacities.
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) PiecewiseLinearMomentCapacityConstitutiveLaw : public ConstitutiveLaw
{
public:
    using BaseType = ConstitutiveLaw;

    KRATOS_CLASS_POINTER_DEFINITION(PiecewiseLinearMomentCapacityConstitutiveLaw);

    // Geometry constants
    static constexpr SizeType strain_size      = 1;
    static constexpr SizeType space_dimenstion = 3;

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    void GetLawFeatures(Features& rFeatures) override;

    double& CalculateValue(ConstitutiveLaw::Parameters& rParameters,
                           const Variable<double>&      rVariable,
                           double&                      rValue) override;
    using BaseType::CalculateValue;

    void CalculateMaterialResponsePK2(Parameters& rValues) override;
    void FinalizeMaterialResponsePK2(Parameters& rValues) override;

    void InitializeMaterial(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector&       rShapeFunctionsValues) override;

    [[nodiscard]] SizeType GetStrainSize() const override;

    [[nodiscard]] int Check(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const ProcessInfo&  rCurrentProcessInfo) const override;

    std::string Info() const override;

private:
    Table<double, double> mStressStrainTable;
    // Unload/reload state (optional, activated when UNRELOAD_MODULUS property is set)
    double mAccumulatedCurvature = 0.0;
    double mPreviousCurvature    = 0.0;
    double mUnReLoadCenter       = 0.0;
    // Optional stored modulus value for unload/reload behavior (set in InitializeMaterial)
    double mUnReLoadModulus = 0.0;

    [[nodiscard]] double CalculateUnReLoadAmplitude() const;
    [[nodiscard]] bool   IsWithinUnReLoading(double Curvature) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos
