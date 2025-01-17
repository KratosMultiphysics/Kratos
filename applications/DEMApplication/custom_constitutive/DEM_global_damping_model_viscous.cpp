//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// Project includes
#include "DEM_global_damping_model_viscous.h"
#include "custom_utilities/GeometryFunctions.h"

namespace Kratos
{
    DEMGlobalDampingModel::Pointer DEMGlobalDampingModelViscous::Clone() const {
        DEMGlobalDampingModel::Pointer p_clone(new DEMGlobalDampingModelViscous(*this));
        return p_clone;
    }

    std::unique_ptr<DEMGlobalDampingModel> DEMGlobalDampingModelViscous::CloneUnique() {
        return Kratos::make_unique<DEMGlobalDampingModelViscous>();
    }

    void DEMGlobalDampingModelViscous::AddGlobalDampingForceAndMoment(SphericParticle* p_element, array_1d<double,3>& total_forces, array_1d<double,3>& total_moment) {
        KRATOS_TRY

        if (p_element->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER) && !p_element->Is(DEMFlags::CUMULATIVE_ZONE)) {
            const array_1d<double, 3>& vel = p_element->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            const double vel_magnitude = DEM_MODULUS_3(vel);
            const double m = p_element->GetMass();
            const double r = p_element->GetRadius();
            const double E = p_element->GetYoung();

            if (vel_magnitude != 0.0) {
                const array_1d<double, 3> damping_force = -2.0 * mGlobalDamping * sqrt(m*r*E) * vel;
                noalias(total_forces) += damping_force;
            }
        }

        KRATOS_CATCH("")
  }
}
