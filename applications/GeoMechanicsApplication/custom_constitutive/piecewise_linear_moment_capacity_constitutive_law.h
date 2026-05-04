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

    void InitializeMaterial(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector&       rShapeFunctionsValues) override;

    [[nodiscard]] SizeType GetStrainSize() const override;

    [[nodiscard]] int Check(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const ProcessInfo&  rCurrentProcessInfo) const override;

    std::string Info() const override;

private:
    [[nodiscard]] double CalculateMomentResponse(double AbsCurvature) const;

    Table<double, double> mStressStrainTable;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos
