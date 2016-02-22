//
// Author: Miguel Angel Celigueta    maceli@cimne.upc.edu
//

// System includes
#include <string>
#include <iostream>
#include <iomanip> // to improve std::cout precision

// Project includes
#include "nanoparticle.h"

namespace Kratos {
void NanoParticle::CustomInitialize()
{
    mSearchRadius = GetGeometry()[0].GetSolutionStepValue(RADIUS);
    mRadius = 1e-7;
    double density            = GetDensity();
    double added_mass_coefficient = 1;
    mRealMass                 = added_mass_coefficient * density * GetVolume();
    GetGeometry()[0].GetSolutionStepValue(NODAL_MASS) = mRealMass;
}

void NanoParticle::ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force,
                             array_1d<double, 3>& additionally_applied_moment,
                             ProcessInfo& r_current_process_info,
                             const array_1d<double,3>& gravity){
    KRATOS_TRY

    array_1d<double, 3> brownian_motion_force; brownian_motion_force.clear();
    array_1d<double, 3> van_der_waals_force; van_der_waals_force.clear();
    array_1d<double, 3> double_layer_force; double_layer_force.clear();

    noalias(additionally_applied_force) += brownian_motion_force + van_der_waals_force + double_layer_force;


    //Now add the contribution of base class function (gravity or other forces added in upper levels):
    SphericParticle::ComputeAdditionalForces(additionally_applied_force, additionally_applied_moment, r_current_process_info, gravity);

    KRATOS_CATCH( "" )
}

double NanoParticle::GetVolume()
{
    return KRATOS_M_PI * mRadius * mRadius * mRadius * mThicknessOverRadius;
}

double NanoParticle::GetSearchRadius()
{
    return mSearchRadius;
}
    
} // namespace Kratos
