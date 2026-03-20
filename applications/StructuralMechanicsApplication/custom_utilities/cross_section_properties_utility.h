// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Siemer
//

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "includes/model_part.h"

namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) CrossSectionPropertiesUtility
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CrossSectionPropertiesUtility);

    CrossSectionPropertiesUtility() = delete;

    static void CalculateBeamProperties(ModelPart& rModelPart);
};

///@}
} // namespace Kratos
