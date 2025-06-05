// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors: Richard Faasse
//

#include "apply_initial_stress_field.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

ApplyInitialUniformStressField::ApplyInitialUniformStressField(ModelPart& rModelPart, const Parameters& rParameters)
    : mrModelPart(rModelPart), mImposedStressVector(rParameters["value"].GetVector())
{
}

void ApplyInitialUniformStressField::ExecuteInitialize()
{
    block_for_each(mrModelPart.Elements(), [this](Element& rElement) {
        std::vector<Vector> dummy_stress_vector(
            rElement.GetGeometry().IntegrationPointsNumber(rElement.GetIntegrationMethod()), mImposedStressVector);
        rElement.SetValuesOnIntegrationPoints(CAUCHY_STRESS_VECTOR, dummy_stress_vector,
                                              mrModelPart.GetProcessInfo());
    });
}

} // namespace Kratos
