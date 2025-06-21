//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// Project includes
#include "DEM_global_damping_nonviscous_constantforcedir.h"
#include "custom_utilities/GeometryFunctions.h"

namespace Kratos
{
    DEMGlobalDampingModel::Pointer DEMGlobalDampingNonViscousCstForceDir::Clone() const {
        DEMGlobalDampingModel::Pointer p_clone(new DEMGlobalDampingNonViscousCstForceDir(*this));
        return p_clone;
    }

    std::unique_ptr<DEMGlobalDampingModel> DEMGlobalDampingNonViscousCstForceDir::CloneUnique() {
        return Kratos::make_unique<DEMGlobalDampingNonViscousCstForceDir>();
    }

    void DEMGlobalDampingNonViscousCstForceDir::AddGlobalDampingForceAndMoment(SphericParticle* p_element, array_1d<double,3>& total_forces, array_1d<double,3>& total_moment) {
        KRATOS_TRY

        const auto& central_node = p_element->GetGeometry()[0];
        const array_1d<double, 3> velocity = central_node.FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3> angular_velocity = central_node.FastGetSolutionStepValue(ANGULAR_VELOCITY);

        if (central_node.IsNot(DEMFlags::FIXED_VEL_X)) {
            total_forces[0] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_forces[0] * velocity[0]));
        }
        if (central_node.IsNot(DEMFlags::FIXED_VEL_Y)) {
            total_forces[1] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_forces[1] * velocity[1]));
        }
        if (central_node.IsNot(DEMFlags::FIXED_VEL_Z)) {
            total_forces[2] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_forces[2] * velocity[2]));
        }

        if (central_node.IsNot(DEMFlags::FIXED_ANG_VEL_X)) {
            total_moment[0] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_moment[0] * angular_velocity[0]));
        }
        if (central_node.IsNot(DEMFlags::FIXED_ANG_VEL_Y)) {
            total_moment[1] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_moment[1] * angular_velocity[1]));
        }
        if (central_node.IsNot(DEMFlags::FIXED_ANG_VEL_Z)) {
            total_moment[2] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_moment[2] * angular_velocity[2]));
        }

        KRATOS_CATCH("")
  }
}
