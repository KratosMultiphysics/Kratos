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

#include "custom_constitutive/linear_elastic_law.h"
#include "geo_mechanics_application_constants.h"

namespace Kratos
{

class ConstitutiveLawDimension;

/**
 * @class GeoLinearElasticPlaneStrain2DLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines a small deformation linear elastic constitutive model for plane strain cases
 * @author Vahid Galavi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoLinearElasticPlaneStrain2DLaw : public GeoLinearElasticLaw
{
public:
    using BaseType = GeoLinearElasticLaw;
    using SizeType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(GeoLinearElasticPlaneStrain2DLaw);
    GeoLinearElasticPlaneStrain2DLaw();

    explicit GeoLinearElasticPlaneStrain2DLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension);
    GeoLinearElasticPlaneStrain2DLaw(const GeoLinearElasticPlaneStrain2DLaw& rOther);
    GeoLinearElasticPlaneStrain2DLaw& operator=(const GeoLinearElasticPlaneStrain2DLaw& rOther);

    GeoLinearElasticPlaneStrain2DLaw(GeoLinearElasticPlaneStrain2DLaw&& rOther)            = delete;
    GeoLinearElasticPlaneStrain2DLaw& operator=(GeoLinearElasticPlaneStrain2DLaw&& rOther) = delete;
    ~GeoLinearElasticPlaneStrain2DLaw() override;

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    bool RequiresInitializeMaterialResponse() override;
    void InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    bool RequiresFinalizeMaterialResponse() override;
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures: The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Dimension of the law:
     * @return The dimension for which the law is working
     */
    SizeType WorkingSpaceDimension() override;

    /**
     * @brief Voigt tensor size:
     * @return The size of the strain vector in Voigt notation
     */
    [[nodiscard]] SizeType GetStrainSize() const override;

    bool IsIncremental() override;

    /**
     * @brief  It returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue) override;
    using ConstitutiveLaw::GetValue;

protected:
    /**
     * @brief It calculates the constitutive matrix C
     * @param C The constitutive matrix
     * @param rValues Parameters of the constitutive law
     */
    void CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief It calculates the stress vector
     * @param rStrainVector The strain vector in Voigt notation
     * @param rStressVector The stress vector in Voigt notation
     * @param rValues Parameters of the constitutive law
     */
    void CalculatePK2Stress(const Vector&                rStrainVector,
                            Vector&                      rStressVector,
                            ConstitutiveLaw::Parameters& rValues) override;

    ///@}

private:
    std::unique_ptr<ConstitutiveLawDimension> mpConstitutiveDimension;
    Vector                                    mStressVector;
    Vector                                    mStressVectorFinalized;
    Vector                                    mDeltaStrainVector;
    Vector                                    mStrainVectorFinalized;
    bool                                      mIsModelInitialized = false;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
}; // Class GeoLinearElasticPlaneStrain2DLaw

} // namespace Kratos