// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/beam_constitutive_law.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

BeamConstitutiveLaw::BeamConstitutiveLaw()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

BeamConstitutiveLaw::BeamConstitutiveLaw(const BeamConstitutiveLaw& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer BeamConstitutiveLaw::Clone() const
{
    BeamConstitutiveLaw::Pointer p_clone(new BeamConstitutiveLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

BeamConstitutiveLaw::~BeamConstitutiveLaw()
{
    // TODO: Add if necessary
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void BeamConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the strain size
    rFeatures.mStrainSize = 1;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}

//************************************************************************************
//************************************************************************************

int BeamConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0) << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;
    const double nu = rMaterialProperties[POISSON_RATIO];
    KRATOS_ERROR_IF((nu > 0.499 && nu < 0.501) || (nu < -0.999 && nu > -1.01)) << "POISSON_RATIO has Key zero or invalid value " << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] <= 0.0) << "DENSITY has Key zero or invalid value " << std::endl;

    return 0;
}

} // Namespace Kratos
