//
// Author: Guillermo Casas gcasas@cimne.upc.edu
//
// Project includes
#include "includes/dem_variables.h"
#include "hybrid_bashforth_scheme.h"

namespace Kratos {

    void HybridBashforthScheme::UpdateTranslationalVariables(
            int StepFlag,
            Node < 3 >& i,
            array_1d<double, 3 >& coor,
            array_1d<double, 3 >& displ,
            array_1d<double, 3 >& delta_displ,
            array_1d<double, 3 >& vel,
            const array_1d<double, 3 >& initial_coor,
            const array_1d<double, 3 >& force,
            const double force_reduction_factor,
            const double mass,
            const double delta_t,
            const bool Fix_vel[3]) {


        if (StepFlag == 1){
            for (int k = 0; k < 3; k++) {
                delta_displ[k] = 0.5 * delta_t * (3 * vel[k] - mOldVelocity[k]);
                displ[k] += delta_displ[k];
                coor[k] = initial_coor[k] + displ[k];
            } // dimensions
        }

        else {
            noalias(mOldVelocity) = vel;

            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {
                    vel[k] += delta_t * force_reduction_factor * force[k] / mass;
                }
            } // dimensions
        }
    }
} //namespace Kratos
