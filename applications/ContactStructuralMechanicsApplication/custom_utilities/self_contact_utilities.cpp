// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//


// System includes
#include <unordered_set>
#include <unordered_map>

// External includes

// Project includes
#include "custom_utilities/self_contact_utilities.h"
#include "contact_structural_mechanics_application_variables.h"

namespace Kratos
{
namespace SelfContactUtilities
{
void ComputeSelfContactPairing(ModelPart& rModelPart)
{
    KRATOS_TRY

    // Iterate over the conditions
    auto& r_conditions_array = rModelPart.Conditions();
    const int num_conditions = static_cast<int>(r_conditions_array.size());
    const auto it_cond_begin = r_conditions_array.begin();

    for(int i = 0; i < num_conditions; ++i) {
        auto it_cond = it_cond_begin + i;
        auto p_indexes_pairs = it_cond->GetValue(INDEX_MAP);

        if (p_indexes_pairs->size() > 0) {
            // TODO
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

} // namespace SelfContactUtilities
} // namespace Kratos
