//
// Author: Guillermo Casas gcasas@cimne.upc.edu
//
#include "includes/dem_variables.h"
#include "hybrid_bashforth_scheme.h"

namespace Kratos {

//    void HybridBashforthScheme::SetIntegrationSchemeInProperties(Properties::Pointer pProp) const {
//        std::cout << "Assigning HybridBashforthScheme to properties " << pProp->Id() << std::endl;
//        pProp->SetValue(DEM_INTEGRATION_SCHEME_POINTER, this->CloneShared());
//    }

    /*void HybridBashforthScheme::AddSpheresVariables(ModelPart & r_model_part, bool TRotationOption){
        DEMIntegrationScheme::AddSpheresVariables(r_model_part, TRotationOption);
    }

    void HybridBashforthScheme::AddClustersVariables(ModelPart & r_model_part, bool TRotationOption){
        DEMIntegrationScheme::AddClustersVariables(r_model_part, TRotationOption);
    }*/

    void HybridBashforthScheme::UpdateTranslationalVariables(
            int StepFlag,
            Node < 3 > & i,
            array_1d<double, 3 >& coor,
            array_1d<double, 3 >& displ,
            array_1d<double, 3 >& delta_displ,
            array_1d<double, 3 >& vel,
            const array_1d<double, 3 >& initial_coor,
            const array_1d<double, 3 >& force,
            const double force_reduction_factor,
            const double mass,
            const double delta_t,
            const bool Fix_vel[3])
    {   
        array_1d<double, 3 >& old_vel = i.FastGetSolutionStepValue(VELOCITY_OLD);

        if (StepFlag == 1){
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {
                    delta_displ[k] = 0.5 * delta_t * (3 * vel[k] - old_vel[k]);
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            } // dimensions
        }

        else {
            noalias(old_vel) = vel;

            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {
                    vel[k] += delta_t * force_reduction_factor * force[k] / mass;
                } else {
                    delta_displ[k] = delta_t * vel[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            } // dimensions
        }
    }
} //namespace Kratos
