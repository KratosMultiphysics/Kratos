#include "custom_constitutive/volumetric_linear_elastic_2D_law.h"

#pragma once 

namespace Kratos 
{

ConstitutiveLaw::Pointer VolumetricLinearElastic2DLaw::Clone() const
{
    VolumetricLinearElastic2DLaw::Pointer p_clone(new VolumetricLinearElastic2DLaw(*this));
    return p_clone; 
}

void VolumetricLinearElastic2DLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize = 3; 
    rFeatures.mSpaceDimension = 2;
}

int VolumetricLinearElastic2DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS)) << "YOUNG MODULUS is not defined in the properties" <<std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO)) << "POISSON RATIO is not defined in the properties" <<std::endl;
    
    KRATOS_ERROR_IF(rMaterialProperties[POISSON_RATIO] < 0.0) << "POISSON RATIO value is lower than 0.0" <<std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[POISSON_RATIO] < 0.5) << "POISSON RATIO value cannot be more than 0.5" <<std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0) << "DENSITY value is lower than 0.0" <<std::endl;
    return 0;
}

void VolumetricLinearElastic2DLaw::CalculateLinearElasticMatrix(Matrix& rConstitutiveMatrix, const Properties& rMaterialProperties)
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double nu = rMaterialProperties[POISSON_RATIO];

    // Computation of Lamé parameters
    const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu  = E / (2.0 * (1.0 + nu));

    //Plane Strain 
    const double beta = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    rConstitutiveMatrix.clear(); 
    rConstitutiveMatrix(0,0) = beta * (1 - nu);
    rConstitutiveMatrix(0,1) = beta * nu;
    rConstitutiveMatrix(1,0) = beta * nu; 
    rConstitutiveMatrix(1,1) = beta * (1 - nu);
    rConstitutiveMatrix(2,2) = beta * (1 - 2 * nu) / 2;

}

void VolumetricLinearElastic2DLaw::CalculateStress(VectorType& rStressVector, VectorType& rStrainVector, Matrix& rConstitutiveMatrix)
{
    if (rStressVector.size() != rStrainVector.size())
        rStressVector.resize(rStrainVector.size(), false);
    noalias(rStressVector) = prod(rConstitutiveMatrix, rStrainVector);
}

void VolumetricLinearElastic2DLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY
    const auto& r_options = rValues.GetOptions();
    auto& rMaterialProperties = rValues.GetMaterialProperties();

    VectorType& strain_vector = rValues.GetStrainVector();
    VectorType& stress_vector = rValues.GetStressVector();

    const SizeType strain_size = GetStrainSize();

    if (r_options.Is(ConstitutiveLaw::COMPUTE_STRESS)){
        if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){
            Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
            this->CalculateLinearElasticMatrix(rConstitutiveMatrix, rMaterialProperties);
            this->CalculateStress(stress_vector, strain_vector, rConstitutiveMatrix);
        } else {
            Matrix ConstitutiveMatrix(strain_size, strain_size);
            noalias(ConstitutiveMatrix) = ZeroMatrix(strain_size, strain_size);

            this->CalculateLinearElasticMatrix(ConstitutiveMatrix, rMaterialProperties);
            this->CalculateStress(stress_vector, strain_vector, ConstitutiveMatrix);
        }
    } else if (r_options.IsNot(ConstitutiveLaw::COMPUTE_STRESS) && r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){
        Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
        this->CalculateLinearElasticMatrix(rConstitutiveMatrix, rMaterialProperties);
    }

    KRATOS_CATCH("")
}

void VolumetricLinearElastic2DLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    VolumetricLinearElastic2DLaw::CalculateMaterialResponseCauchy(rValues);
}

}