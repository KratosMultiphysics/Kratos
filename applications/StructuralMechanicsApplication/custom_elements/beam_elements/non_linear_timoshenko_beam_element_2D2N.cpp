// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes

// Application includes
#include "non_linear_timoshenko_beam_element_2D2N.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LinearTimoshenkoBeamElement2D2N);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LinearTimoshenkoBeamElement2D2N);
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
