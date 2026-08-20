// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

#include "constitutive_law_dimension.h"
#include "custom_constitutive/linear_elastic_law.h"
#include "custom_constitutive/youngs_modulus_formulations.h"
#include "includes/properties.h"

namespace Kratos
{

/**
 * @class GeoIncrementalLinearElasticLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines an incremental linear elastic constitutive model for plane strain and 3D cases
 * @author Vahid Galavi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoIncrementalLinearElasticLaw : public GeoLinearElasticLaw
{
public:
    using BaseType = GeoLinearElasticLaw;
    using SizeType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(GeoIncrementalLinearElasticLaw);
    GeoIncrementalLinearElasticLaw() = default;

    explicit GeoIncrementalLinearElasticLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension);
    GeoIncrementalLinearElasticLaw(const GeoIncrementalLinearElasticLaw& rOther);
    GeoIncrementalLinearElasticLaw& operator=(const GeoIncrementalLinearElasticLaw& rOther);

    GeoIncrementalLinearElasticLaw(GeoIncrementalLinearElasticLaw&& rOther) noexcept = default;
    GeoIncrementalLinearElasticLaw& operator=(GeoIncrementalLinearElasticLaw&& rOther) noexcept = default;
    ~GeoIncrementalLinearElasticLaw() override = default;

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    bool RequiresInitializeMaterialResponse() override;
    void InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters) override;

    bool RequiresFinalizeMaterialResponse() override;
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters) override;
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rParameters) override;

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures: The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Dimension of the law:
     * @return The dimension for which the law is working
     */
    [[nodiscard]] SizeType WorkingSpaceDimension() override;

    /**
     * @brief Voigt tensor size:
     * @return The size of the strain vector in Voigt notation
     */
    [[nodiscard]] SizeType GetStrainSize() const override;

    bool IsIncremental() override;

    /**
     * @brief  It returns the value of a specified variable
     * @param rVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return The value of the specified variable
     */
    bool& GetValue(const Variable<bool>& rVariable, bool& rValue) override;
    using ConstitutiveLaw::GetValue;

    int Check(const Properties&   rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo&  rCurrentProcessInfo) const override;

    /**
     * @brief It resets all the member variables and flags
     */
    void ResetMaterial(const Properties&, const GeometryType&, const Vector&) override;

protected:
    /**
     * @brief It calculates the constitutive matrix
     * @param rElasticMatrix The constitutive matrix
     * @param rParameters Parameters of the constitutive law
     */
    void CalculateElasticMatrix(Matrix& rElasticMatrix, ConstitutiveLaw::Parameters& rParameters) override;

    /**
     * @brief It calculates the stress vector
     * @param rStrainVector The strain vector in Voigt notation
     * @param rStressVector The stress vector in Voigt notation
     * @param rParameters Parameters of the constitutive law
     */
    void CalculatePK2Stress(const Vector&                rStrainVector,
                            Vector&                      rStressVector,
                            ConstitutiveLaw::Parameters& rParameters) override;

    ///@}

private:
    std::unique_ptr<ConstitutiveLawDimension> mpConstitutiveDimension;
    Vector                                    mStressVector;
    Vector                                    mStressVectorFinalized;
    Vector                                    mDeltaStrainVector;
    Vector                                    mStrainVectorFinalized;
    bool                                      mIsModelInitialized = false;
    Formulations::YoungsModulusVariant        mFormulation;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
}; // Class GeoIncrementalLinearElasticLaw

} // namespace Kratos
