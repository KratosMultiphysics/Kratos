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
#include "includes/properties.h"
#include "custom_constitutive/truss_constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TrussConstitutiveLaw::TrussConstitutiveLaw()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

TrussConstitutiveLaw::TrussConstitutiveLaw(const TrussConstitutiveLaw& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer TrussConstitutiveLaw::Clone() const
{
    return Kratos::make_shared<TrussConstitutiveLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

TrussConstitutiveLaw::~TrussConstitutiveLaw()
{
    // TODO: Add if necessary
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void TrussConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the strain size
    rFeatures.mStrainSize = 1;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}
//************************************************************************************
//************************************************************************************

array_1d<double, 3 > & TrussConstitutiveLaw::GetValue(
    const Variable<array_1d<double, 3 > >& rThisVariable,
    array_1d<double, 3 > & rValue)
{
    KRATOS_ERROR << "Can't get the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

double& TrussConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    if(rThisVariable == TANGENT_MODULUS) rValue = rParameterValues.GetMaterialProperties()[YOUNG_MODULUS];
    else KRATOS_ERROR << "Can't calculate the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

Vector& TrussConstitutiveLaw::CalculateValue(
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

array_1d<double, 3 > & TrussConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<array_1d<double, 3 > >& rVariable,
	array_1d<double, 3 > & rValue)
    {
        if (rVariable == FORCE)
        {
            constexpr SizeType dimension = 3;
            rValue = ZeroVector(dimension);
            //rValue[0] = this->mStressState;
            rValue[0] = this->CalculateStressElastic(rParameterValues);
            rValue[1] = 0.0;
            rValue[2] = 0.0;
        }
        else KRATOS_ERROR << "Can't calculate the specified value" << std::endl;
        return rValue;
    }

//************************************************************************************
//************************************************************************************

void TrussConstitutiveLaw::CalculateMaterialResponse(
    const Vector& rStrainVector,const Matrix& rDeformationGradient,
    Vector& rStressVector,Matrix& rAlgorithmicTangent,
    const ProcessInfo& rCurrentProcessInfo,const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues,
    bool CalculateStresses,int CalculateTangent,bool SaveInternalVariables)
{
    const double axial_strain = rStrainVector[0];
    const double youngs_modulus = rMaterialProperties[YOUNG_MODULUS];

    if (rStressVector.size() != 1) rStressVector.resize(1);
    rStressVector[0] = youngs_modulus*axial_strain;
}

//************************************************************************************
//************************************************************************************

double TrussConstitutiveLaw::CalculateStressElastic(
    ConstitutiveLaw::Parameters& rParameterValues) const
{
    Vector current_strain = ZeroVector(1);
    rParameterValues.GetStrainVector(current_strain);

    const double current_stress =
     rParameterValues.GetMaterialProperties()[YOUNG_MODULUS]*current_strain[0];
    return current_stress;
}

//************************************************************************************
//************************************************************************************

int TrussConstitutiveLaw::Check(
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
