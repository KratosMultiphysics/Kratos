// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder
//

// System includes

// External includes

// Project includes
#include "generalized_influence_functions_process.h"
#include "custom_response_functions/adjoint_elements/generalized_influence_functions_extension.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

void GeneralizedInfluenceFunctionsProcess::ExecuteFinalizeSolutionStep()
{
    for(int i=0; i< static_cast<int>(mrModelPart.Elements().size()); ++i)
    {
        auto it = mrModelPart.ElementsBegin() + i;
        it->SetValue(INFLUENCE_FUNCTIONS_EXTENSIONS, Kratos::make_shared<GeneralizedInfluenceFunctionsExtension>(mSettings));
    }
}

}  // namespace Kratos.



