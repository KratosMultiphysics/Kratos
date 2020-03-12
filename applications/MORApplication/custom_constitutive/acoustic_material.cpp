// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: 
//
//  Main authors:    
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/acoustic_material.h"
#include "includes/checks.h"

#include "mor_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

AcousticMaterial::AcousticMaterial()
    : ConstitutiveLaw()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

AcousticMaterial::AcousticMaterial(const AcousticMaterial& rOther)
    : ConstitutiveLaw(rOther)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer AcousticMaterial::Clone() const
{
    return Kratos::make_shared<AcousticMaterial>(*this);
}

/*******************************DESTRUCTOR******************************************/
/***********************************************************************************/

AcousticMaterial::~AcousticMaterial()
{
};


/***********************************************************************************/
/***********************************************************************************/

void AcousticMaterial::CheckClearElasticMatrix(Matrix& rConstitutiveMatrix)
{
    const SizeType size_system = this->GetStrainSize();
    if (rConstitutiveMatrix.size1() != size_system || rConstitutiveMatrix.size2() != size_system)
        rConstitutiveMatrix.resize(size_system, size_system, false);
    rConstitutiveMatrix.clear();
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticMaterial::CalculateElasticMatrix(
    Matrix& rConstitutiveMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double RO = r_material_properties[DENSITY];

    this->CheckClearElasticMatrix(rConstitutiveMatrix);

    const double c1 = sqrt(E/RO);

    rConstitutiveMatrix( 0, 0 ) = c1;

   
}

/***********************************************************************************/
/***********************************************************************************/


} // Namespace Kratos
