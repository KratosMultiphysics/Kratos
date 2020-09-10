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
#include "custom_advanced_constitutive/hyper_elastic_isotropic_henky_1d.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticIsotropicHenky1D::HyperElasticIsotropicHenky1D()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticIsotropicHenky1D::HyperElasticIsotropicHenky1D(const HyperElasticIsotropicHenky1D& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticIsotropicHenky1D::Clone() const
{
    return Kratos::make_shared<HyperElasticIsotropicHenky1D>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticIsotropicHenky1D::~HyperElasticIsotropicHenky1D()
{
    // TODO: Add if necessary
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HyperElasticIsotropicHenky1D::GetLawFeatures(Features& rFeatures)
{
    //Set the strain size
    rFeatures.mStrainSize = 1;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}
//************************************************************************************
//************************************************************************************

array_1d<double, 3 > & HyperElasticIsotropicHenky1D::GetValue(
    const Variable<array_1d<double, 3 > >& rThisVariable,
    array_1d<double, 3 > & rValue)
{
    KRATOS_ERROR << "Can't get the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

double& HyperElasticIsotropicHenky1D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    if(rThisVariable == TANGENT_MODULUS)
    {
        const double E(rParameterValues.GetMaterialProperties()[YOUNG_MODULUS]);

        Vector current_strain = ZeroVector(1);
        rParameterValues.GetStrainVector(current_strain);
        const double E_11(current_strain[0]);

        rValue = (E-(E*std::log(2.0*E_11 + 1.0)))/std::pow(2.0*E_11 + 1.0, 2);
;
    }
    else KRATOS_ERROR << "Can't calculate the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

Vector& HyperElasticIsotropicHenky1D::CalculateValue(
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

array_1d<double, 3 > & HyperElasticIsotropicHenky1D::CalculateValue(
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
void HyperElasticIsotropicHenky1D::CalculateMaterialResponsePK2(Parameters& rValues)
{
    Vector& stress_vector = rValues.GetStressVector();
    if (stress_vector.size() != 1) stress_vector.resize(1, false);
    stress_vector[0] = this->CalculateStressElastic(rValues);
}
//************************************************************************************
//************************************************************************************

double HyperElasticIsotropicHenky1D::CalculateStressElastic(
    ConstitutiveLaw::Parameters& rParameterValues) const
{
    const double E(rParameterValues.GetMaterialProperties()[YOUNG_MODULUS]);

    Vector current_strain = ZeroVector(1);
    rParameterValues.GetStrainVector(current_strain);
    const double E_11(current_strain[0]);

    const double current_stress = 1.0*E*std::log(2.0*E_11 + 1.0)/(4.0*E_11 + 2.0);

    return current_stress;
}

//************************************************************************************
//************************************************************************************

int HyperElasticIsotropicHenky1D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_ERROR_IF(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS] <= 0.00)
     << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;

    KRATOS_ERROR_IF(DENSITY.Key() == 0 || rMaterialProperties[DENSITY] < 0.00)
     << "DENSITY has Key zero or invalid value " << std::endl;

    return 0;

}

} // Namespace Kratos
