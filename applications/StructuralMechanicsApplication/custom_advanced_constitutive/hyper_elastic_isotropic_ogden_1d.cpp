// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter
//
// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/properties.h"
#include "custom_advanced_constitutive/hyper_elastic_isotropic_ogden_1d.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticIsotropicOgden1D::HyperElasticIsotropicOgden1D()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticIsotropicOgden1D::HyperElasticIsotropicOgden1D(const HyperElasticIsotropicOgden1D& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticIsotropicOgden1D::Clone() const
{
    return Kratos::make_shared<HyperElasticIsotropicOgden1D>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticIsotropicOgden1D::~HyperElasticIsotropicOgden1D()
{
    // TODO: Add if necessary
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HyperElasticIsotropicOgden1D::GetLawFeatures(Features& rFeatures)
{
    //Set the strain size
    rFeatures.mStrainSize = 1;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}
//************************************************************************************
//************************************************************************************

array_1d<double, 3 > & HyperElasticIsotropicOgden1D::GetValue(
    const Variable<array_1d<double, 3 > >& rThisVariable,
    array_1d<double, 3 > & rValue)
{
    KRATOS_ERROR << "Can't get the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

double& HyperElasticIsotropicOgden1D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    if(rThisVariable == TANGENT_MODULUS)
    {
        const double E(rParameterValues.GetMaterialProperties()[YOUNG_MODULUS]);
        const double beta_1(rParameterValues.GetMaterialProperties()[OGDEN_BETA_1]);
        const double beta_2(rParameterValues.GetMaterialProperties()[OGDEN_BETA_2]);

        Vector current_strain = ZeroVector(1);
        rParameterValues.GetStrainVector(current_strain);
        const double E_11(current_strain[0]);

        rValue = E*(1.0*beta_1*std::pow(2.0*E_11 + 1.0, (1.0/2.0)*beta_1)/std::pow(2.0*E_11 + 1.0, 2) -
         1.0*beta_2*std::pow(2.0*E_11 + 1.0, (1.0/2.0)*beta_2)/std::pow(2.0*E_11 + 1.0, 2) -
         2.0*std::pow(2.0*E_11 + 1.0, (1.0/2.0)*beta_1)/std::pow(2.0*E_11 + 1.0, 2) +
         2.0*std::pow(2.0*E_11 + 1.0, (1.0/2.0)*beta_2)/std::pow(2.0*E_11 + 1.0, 2))/(beta_1 - beta_2);
    }
    else KRATOS_ERROR << "Can't calculate the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

Vector& HyperElasticIsotropicOgden1D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    if(rThisVariable == NORMAL_STRESS)
    {
        const double current_stress = this->CalculateStressElastic(rParameterValues);
        constexpr SizeType dofs = 6;
        rValue = ZeroVector(dofs);
        rValue[0] = -1.0 * current_stress;
        rValue[3] = 1.0 * current_stress;
    }
    else KRATOS_ERROR << "Can't calculate the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

array_1d<double, 3 > & HyperElasticIsotropicOgden1D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<array_1d<double, 3 > >& rVariable,
	array_1d<double, 3 > & rValue)
    {
        if (rVariable == FORCE)
        {
            constexpr SizeType dimension = 3;
            rValue = ZeroVector(dimension);
            rValue[0] = this->CalculateStressElastic(rParameterValues);
            rValue[1] = 0.0;
            rValue[2] = 0.0;
        }
        else KRATOS_ERROR << "Can't calculate the specified value" << std::endl;
        return rValue;
    }

//************************************************************************************
//************************************************************************************
void HyperElasticIsotropicOgden1D::CalculateMaterialResponsePK2(Parameters& rValues)
{
    Vector& stress_vector = rValues.GetStressVector();
    if (stress_vector.size() != 1) stress_vector.resize(1, false);
    stress_vector[0] = this->CalculateStressElastic(rValues);
}
//************************************************************************************
//************************************************************************************

double HyperElasticIsotropicOgden1D::CalculateStressElastic(
    ConstitutiveLaw::Parameters& rParameterValues) const
{
    const double E(rParameterValues.GetMaterialProperties()[YOUNG_MODULUS]);
    const double beta_1(rParameterValues.GetMaterialProperties()[OGDEN_BETA_1]);
    const double beta_2(rParameterValues.GetMaterialProperties()[OGDEN_BETA_2]);

    Vector current_strain = ZeroVector(1);
    rParameterValues.GetStrainVector(current_strain);
    const double E_11(current_strain[0]);

    const double current_stress = E*(1.0*std::pow(2.0*E_11 + 1.0, (1.0/2.0)*beta_1)/(2.0*E_11 + 1.0) -
     1.0*std::pow(2.0*E_11 + 1.0, (1.0/2.0)*beta_2)/(2.0*E_11 + 1.0))/(beta_1 - beta_2);

    return current_stress;
}

//************************************************************************************
//************************************************************************************

int HyperElasticIsotropicOgden1D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_ERROR_IF(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS] <= 0.00)
     << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(OGDEN_BETA_1);
    KRATOS_CHECK(rMaterialProperties.Has(OGDEN_BETA_1));

    KRATOS_CHECK_VARIABLE_KEY(OGDEN_BETA_2);
    KRATOS_CHECK(rMaterialProperties.Has(OGDEN_BETA_2));

    KRATOS_ERROR_IF(rMaterialProperties[OGDEN_BETA_1]==rMaterialProperties[OGDEN_BETA_2]) << "ogden parameters must not be the same" << std::endl;

    KRATOS_ERROR_IF(DENSITY.Key() == 0 || rMaterialProperties[DENSITY] < 0.00)
     << "DENSITY has Key zero or invalid value " << std::endl;

    return 0;

}

} // Namespace Kratos
