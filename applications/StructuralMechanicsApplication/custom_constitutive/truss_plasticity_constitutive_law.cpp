// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/truss_plasticity_constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TrussPlasticityConstitutiveLaw::TrussPlasticityConstitutiveLaw()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

TrussPlasticityConstitutiveLaw::TrussPlasticityConstitutiveLaw(const TrussPlasticityConstitutiveLaw& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer TrussPlasticityConstitutiveLaw::Clone() const
{
    TrussPlasticityConstitutiveLaw::Pointer p_clone(new TrussPlasticityConstitutiveLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

TrussPlasticityConstitutiveLaw::~TrussPlasticityConstitutiveLaw()
{
    // TODO: Add if necessary
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void TrussPlasticityConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the strain size
    rFeatures.mStrainSize = 1;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}

//************************************************************************************
//************************************************************************************

int TrussPlasticityConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));

    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK(rMaterialProperties.Has(DENSITY));

    KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS);
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS));

    KRATOS_CHECK_VARIABLE_KEY(HARDENING_MODULUS_1D);
    KRATOS_CHECK(rMaterialProperties.Has(HARDENING_MODULUS_1D));
    return 0;
}

} // Namespace Kratos
