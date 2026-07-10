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

#include "apply_initial_uniform_stress_field.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
using namespace std::string_literals;

ApplyInitialUniformStressField::ApplyInitialUniformStressField(ModelPart& rModelPart, const Parameters& rParameters)
    : mrModelPart(rModelPart), mImposedStressVector(rParameters["value"].GetVector())
{
    block_for_each(mrModelPart.Elements(), [this](const auto& rElement) {
        const auto p_constitutive_law = rElement.GetProperties()[CONSTITUTIVE_LAW];
        KRATOS_ERROR_IF_NOT(mImposedStressVector.size() == p_constitutive_law->GetStrainSize())
            << "The size of the input stress vector for applying a uniform initial "
               "stress field must match the strain size of the constitutive law, which is "
            << p_constitutive_law->GetStrainSize() << ", but is " << mImposedStressVector.size()
            << " for element " << rElement.Id() << " in model part '" << mrModelPart.Name()
            << "'. Please check the process parameters.\n";
    });
}

void ApplyInitialUniformStressField::ExecuteInitialize()
{
    block_for_each(mrModelPart.Elements(), [this](auto& rElement) {
        std::vector<Vector> imposed_stress_vectors(
            rElement.GetGeometry().IntegrationPointsNumber(rElement.GetIntegrationMethod()), mImposedStressVector);
        rElement.SetValuesOnIntegrationPoints(CAUCHY_STRESS_VECTOR, imposed_stress_vectors,
                                              mrModelPart.GetProcessInfo());
    });
}

std::string ApplyInitialUniformStressField::Info() const
{
    return "ApplyInitialUniformStressField"s;
}

} // namespace Kratos
