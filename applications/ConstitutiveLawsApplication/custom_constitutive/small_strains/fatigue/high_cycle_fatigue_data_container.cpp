// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/small_strains/fatigue/high_cycle_fatigue_data_container.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

HighCycleFatigueDataContainer::HighCycleFatigueDataContainer()
    : ElasticIsotropic3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

HighCycleFatigueDataContainer::HighCycleFatigueDataContainer(const HighCycleFatigueDataContainer& rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer HighCycleFatigueDataContainer::Clone() const
{
    return Kratos::make_shared<HighCycleFatigueDataContainer>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

HighCycleFatigueDataContainer::~HighCycleFatigueDataContainer()
{
}

} // Namespace Kratos
