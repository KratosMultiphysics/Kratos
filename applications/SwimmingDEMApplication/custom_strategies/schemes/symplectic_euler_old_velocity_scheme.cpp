//
// Author: Guillermo Casas gcasas@cimne.upc.edu
//
#include "includes/dem_variables.h"
#include "symplectic_euler_old_velocity_scheme.h"

namespace Kratos {

void SymplecticEulerOldVelocityScheme::UpdateTranslationalVariables(
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
        const bool Fix_vel[3]) {

    array_1d<double, 3 >& old_vel = i.FastGetSolutionStepValue(VELOCITY_OLD);
    noalias(old_vel) = vel;

    for (int k = 0; k < 3; k++) {
        if (Fix_vel[k] == false) {
            vel[k] += delta_t * force_reduction_factor * force[k] / mass;
            delta_displ[k] = delta_t * vel[k];
            displ[k] += delta_displ[k];
            coor[k] = initial_coor[k] + displ[k];
        } else {
            delta_displ[k] = delta_t * vel[k];
            displ[k] += delta_displ[k];
            coor[k] = initial_coor[k] + displ[k];
        }
    } // dimensions
}

} //namespace Kratos
