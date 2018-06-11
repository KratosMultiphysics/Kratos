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
#include <iostream>

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
/*     TrussConstitutiveLaw::Pointer p_clone(new TrussConstitutiveLaw(*this));
    return p_clone;
 */
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

double& TrussConstitutiveLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if(rThisVariable == VON_MISES_STRESS_MIDDLE_SURFACE) rValue = this->mStressState;
    else KRATOS_ERROR << "can't get the specified value" << std::endl;
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
    else KRATOS_ERROR << "can't calculate the specified value" << std::endl;
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
        const SizeType dofs = 6;
        if(rValue.size()!=dofs) rValue.resize(dofs);
        rValue[0] = -1.0 * this->mStressState;
        rValue[3] = 1.0 * this->mStressState;
    }
    else KRATOS_ERROR << "can't calculate the specified value" << std::endl;
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

    //if (SaveInternalVariables) this->mStressState = rStressVector[0];
    this->mStressState = rStressVector[0];
}


//************************************************************************************
//************************************************************************************

int TrussConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
)
{
    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS] <= 0.00)
    {
        KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;
    }

    if(DENSITY.Key() == 0 || rMaterialProperties[DENSITY] < 0.00)
    {
        KRATOS_ERROR << "DENSITY has Key zero or invalid value " << std::endl;
    }

    return 0;

}

} // Namespace Kratos
