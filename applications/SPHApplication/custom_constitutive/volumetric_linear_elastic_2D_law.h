#pragma once 

#include "includes/constitutive_law.h"
#include "sph_application_variables.h"

namespace Kratos
{

/**
 * @class VolumetricLinearElasticLaw
 * @ingroup SPHApplication
 * @brief This class is used by the SPH particle element
 * @details in a second moment constitutive laws from CLapp could be used directly
 * @author
 */

class KRATOS_API(SPH_APPLICATION) VolumetricLinearElastic2DLaw
    : public ConstitutiveLaw
{

public:

    using SizeType = std::size_t;
    using VectorType = Vector;

    KRATOS_CLASS_POINTER_DEFINITION(VolumetricLinearElastic2DLaw);

    VolumetricLinearElastic2DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override;

    VolumetricLinearElastic2DLaw(const VolumetricLinearElastic2DLaw &other)
    {
    }

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeature
     */
    void GetLawFeatures(Features& rFeatures) override;

    SizeType GetStrainSize() const override
    {
        return 3;
    }

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * @param rMaterialProperties: The properties of the material
     * @param rElementGeometry: The geometry of the element
     * @param rCurrentProcessInfo: The current process info instance
     */
    int Check(const Properties& rMaterialProperties, 
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    /**
     * @brief
     */
    virtual void CalculateStress(VectorType& rStressVector, VectorType& rStrainVector, Matrix& rConstitutiveMatrix);
    
    /**
     * @brief Computes the Linear Elastic Constitutive matrix 
     */
    virtual void CalculateLinearElasticMatrix(Matrix& rConstitutiveMatrix, const Properties& rMaterialProperties);

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;
    
    /**
     * @brief Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues) override;

protected:

private:

};

}