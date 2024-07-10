// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//
#include "reset_displacement_process.h"
#include "includes/initial_state.h"
#include "includes/mat_variables.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"

#include <boost/format/free_funcs.hpp>

namespace Kratos
{

ResetDisplacementProcess::ResetDisplacementProcess(ModelPart& rModelPart, const Parameters&)
    : mrModelPart(rModelPart)
{
}

void ResetDisplacementProcess::ExecuteInitialize()
{
    block_for_each(mrModelPart.Elements(), [this](Element& rElement) {
        std::vector<Vector> stresses_on_integration_points;
        rElement.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stresses_on_integration_points,
                                              mrModelPart.GetProcessInfo());

        std::vector<ConstitutiveLaw::Pointer> constitutive_laws;
        rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws, mrModelPart.GetProcessInfo());

        for (auto i = std::size_t{0}; i < constitutive_laws.size(); ++i) {
            auto p_initial_state = make_intrusive<InitialState>();
            p_initial_state->SetInitialStressVector(stresses_on_integration_points[i]);
            constitutive_laws[i]->SetInitialState(p_initial_state);
        }
    });
}

} // namespace Kratos