// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alireza Taherzadeh Fard
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/cohesive_law.h"

#include "constitutive_laws_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

CohesiveLaw::CohesiveLaw()
    : ElasticIsotropic3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

CohesiveLaw::CohesiveLaw(const CohesiveLaw& rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer CohesiveLaw::Clone() const
{
    return Kratos::make_shared<CohesiveLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

CohesiveLaw::~CohesiveLaw()
{
}


} // Namespace Kratos
