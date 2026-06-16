// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#pragma once

#include "constitutive_law_dimension.h"
#include "custom_constitutive/linear_elastic_law.h"

namespace Kratos
{
/**
 * @class GeoIncrementalLinearElasticEurLaw
 * @ingroup GeoMechanicsApplication
 * @brief Incremental linear elastic law with a stress-dependent Young's modulus following an E_ur type formulation
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoIncrementalLinearElasticEurLaw : public GeoLinearElasticLaw
{
public:
    using BaseType = GeoLinearElasticLaw;
    using SizeType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(GeoIncrementalLinearElasticEurLaw);
    GeoIncrementalLinearElasticEurLaw() = default;

    explicit GeoIncrementalLinearElasticEurLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension);
    GeoIncrementalLinearElasticEurLaw(const GeoIncrementalLinearElasticEurLaw& rOther);
    GeoIncrementalLinearElasticEurLaw& operator=(const GeoIncrementalLinearElasticEurLaw& rOther);

    GeoIncrementalLinearElasticEurLaw(GeoIncrementalLinearElasticEurLaw&& rOther) noexcept = default;
    GeoIncrementalLinearElasticEurLaw& operator=(GeoIncrementalLinearElasticEurLaw&& rOther) noexcept = default;
    ~GeoIncrementalLinearElasticEurLaw() override = default;

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    bool RequiresInitializeMaterialResponse() override;
    void InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    bool RequiresFinalizeMaterialResponse() override;
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    void GetLawFeatures(Features& rFeatures) override;

    [[nodiscard]] SizeType WorkingSpaceDimension() override;
    [[nodiscard]] SizeType GetStrainSize() const override;

    bool IsIncremental() override;

    bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue) override;
    using ConstitutiveLaw::GetValue;

    int Check(const Properties&   rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo&  rCurrentProcessInfo) const override;

    void ResetMaterial(const Properties&, const GeometryType&, const Vector&) override;

protected:
    void CalculateElasticMatrix(Matrix& rElasticMatrix, ConstitutiveLaw::Parameters& rValues) override;

    void CalculatePK2Stress(const Vector&                rStrainVector,
                            Vector&                      rStressVector,
                            ConstitutiveLaw::Parameters& rValues) override;

private:
    [[nodiscard]] double CalculateMinorPrincipalEffectiveStress() const;
    [[nodiscard]] double CalculateStressDependentYoungsModulus(const Properties& rProperties) const;

    std::unique_ptr<ConstitutiveLawDimension> mpConstitutiveDimension;
    Vector                                    mStressVector;
    Vector                                    mStressVectorFinalized;
    Vector                                    mDeltaStrainVector;
    Vector                                    mStrainVectorFinalized;
    bool                                      mIsModelInitialized = false;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos
