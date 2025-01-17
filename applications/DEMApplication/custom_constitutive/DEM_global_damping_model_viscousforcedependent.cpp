//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// Project includes
#include "DEM_global_damping_model_viscousforcedependent.h"
#include "custom_utilities/GeometryFunctions.h"

namespace Kratos
{
    DEMGlobalDampingModel::Pointer DEMGlobalDampingModelViscousForceDependent::Clone() const {
        DEMGlobalDampingModel::Pointer p_clone(new DEMGlobalDampingModelViscousForceDependent(*this));
        return p_clone;
    }

    std::unique_ptr<DEMGlobalDampingModel> DEMGlobalDampingModelViscousForceDependent::CloneUnique() {
        return Kratos::make_unique<DEMGlobalDampingModelViscousForceDependent>();
    }

    void DEMGlobalDampingModelViscousForceDependent::AddGlobalDampingForceAndMoment(SphericParticle* p_element, array_1d<double,3>& total_forces, array_1d<double,3>& total_moment) {
        KRATOS_TRY
        
        const auto& central_node = p_element->GetGeometry()[0];
        array_1d<double, 3> velocity = central_node.FastGetSolutionStepValue(VELOCITY);
        GeometryFunctions::normalize(velocity);
        const double force = DEM_MODULUS_3(total_forces);
        
        if (central_node.IsNot(DEMFlags::FIXED_VEL_X)) {
            total_forces[0] -= mGlobalDamping * force * velocity[0];
        }
        if (central_node.IsNot(DEMFlags::FIXED_VEL_Y)) {
            total_forces[1] -= mGlobalDamping * force * velocity[1];
        }
        if (central_node.IsNot(DEMFlags::FIXED_VEL_Z)) {
            total_forces[2] -= mGlobalDamping * force * velocity[2];
        }

        KRATOS_CATCH("")
  }
}
