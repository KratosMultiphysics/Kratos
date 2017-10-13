// Author: Guillermo Casas, gcasas@cimne.upc.edu

#include "adams_bashforth_strategy.h"

namespace Kratos {

void AdamsBashforthStrategy::ReconstructForces(ModelPart& r_model_part)
{
    KRATOS_TRY

    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    ElementsArrayType& pElements = GetElements(r_model_part);

    if (pElements.size()){
        ElementsArrayType::iterator it_0 = pElements.ptr_begin(); //first element (any element will do)
        ModelPart::NodeType& p_any_node   = it_0->GetGeometry()[0];
        const bool has_virtual_mass_force = p_any_node.SolutionStepsDataHas(VIRTUAL_MASS_FORCE);
        const bool has_basset_force       = p_any_node.SolutionStepsDataHas(BASSET_FORCE);

        if (has_virtual_mass_force || has_basset_force){ // it is worth it to loop

            #pragma omp parallel for
            for (int k = 0; k < (int) pElements.size(); k++){
                ElementsArrayType::iterator it = pElements.ptr_begin() + k;
                ModelPart::NodeType& p_node = it->GetGeometry()[0];

                if (has_virtual_mass_force){
                    array_1d<double, 3>& virtual_mass_force = p_node.FastGetSolutionStepValue(VIRTUAL_MASS_FORCE);
                    it->Calculate(VIRTUAL_MASS_FORCE, virtual_mass_force, r_process_info);
                }

                if (has_basset_force){
                    array_1d<double, 3>& has_basset_force = p_node.FastGetSolutionStepValue(BASSET_FORCE);
                    it->Calculate(BASSET_FORCE, has_basset_force, r_process_info);
                }
            }
        }
    }

    KRATOS_CATCH("")
}

double AdamsBashforthStrategy::Solve() {
    KRATOS_TRY
    ModelPart& r_model_part = GetModelPart();

    if (mFirstStep){
        mFirstStep = false;
        PerformTimeIntegrationOfMotion(1);
    }

    else {
        mFirstStep = true;
        InitializeSolutionStep();
        SearchDEMOperations(r_model_part);
        SearchFEMOperations(r_model_part);
        ForceOperations(r_model_part);
        PerformTimeIntegrationOfMotion(2);
        ReconstructForces(r_model_part);
        FinalizeSolutionStep();
    }
    return 0.00;

    KRATOS_CATCH("")
}

} // namespace Kratos
