//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Chengshun Shang (cshang@cimne.upc.edu)
//                 Joaquin Irazabal (jirazabal@cimne.upc.edu)
//

// Description:
// This model can be found as BRFM in [Joaquin Irazabal, 2017, Numerical modelling of granular materials -
// with spherical discrete particles and the bounded rolling friction model. Application to railway ballast]

#include "DEM_rolling_friction_model_bounded.h"
#include "custom_utilities/GeometryFunctions.h"

namespace Kratos{

    DEMRollingFrictionModel::Pointer DEMRollingFrictionModelBounded::Clone() const {
        DEMRollingFrictionModel::Pointer p_clone(new DEMRollingFrictionModelBounded(*this));
        return p_clone;
    }

    std::unique_ptr<DEMRollingFrictionModel> DEMRollingFrictionModelBounded::CloneUnique() {
        return Kratos::make_unique<DEMRollingFrictionModelBounded>();
    }

    void DEMRollingFrictionModelBounded::Check(Properties::Pointer pProp) const {
        
        if(!pProp->Has(DEM_ROLLING_FRICTION_MODEL_NAME)) 
        {
            KRATOS_ERROR << "Variable DEM_ROLLING_FRICTION_MODEL_NAME should be presented."<<std::endl;
        }
        
    }

    bool DEMRollingFrictionModelBounded::CheckIfThisModelRequiresRecloningForEachNeighbour() {
        
        KRATOS_TRY

        return false;

        KRATOS_CATCH("")
    }

    void DEMRollingFrictionModelBounded::InitializeSolutionStep() {
        
        KRATOS_TRY

        mRollingResistance = 0.0;

        KRATOS_CATCH("")
    }

    void DEMRollingFrictionModelBounded::ComputeRollingResistance(SphericParticle* p_element, SphericParticle* p_neighbor, double LocalContactForce[3]) {

        KRATOS_TRY

        Properties& r_properties = p_element->GetProperties().GetSubProperties(p_neighbor->GetProperties().Id());
        const double min_radius = std::min(p_element->GetRadius(), p_neighbor->GetRadius());
        const double equiv_rolling_friction_coeff = r_properties[ROLLING_FRICTION] * min_radius;

        mRollingResistance += std::abs(LocalContactForce[2]) * equiv_rolling_friction_coeff;

        KRATOS_CATCH("")
    }

    void DEMRollingFrictionModelBounded::ComputeRollingResistanceWithWall(SphericParticle* p_element, Condition* const wall, double LocalContactForce[3]) {

        KRATOS_TRY

        Properties& r_properties = p_element->GetProperties().GetSubProperties(wall->GetProperties().Id());
        const double equiv_rolling_friction_coeff = r_properties[ROLLING_FRICTION] * p_element->GetRadius();

        mRollingResistance += std::abs(LocalContactForce[2]) * equiv_rolling_friction_coeff;

        KRATOS_CATCH("")
    }

    void DEMRollingFrictionModelBounded::DoFinalOperations(SphericParticle* p_element, double dt, array_1d<double, 3>& mContactMoment)
    {
        KRATOS_TRY

        array_1d<double, 3> & rolling_resistance_moment = p_element->GetGeometry()[0].FastGetSolutionStepValue(ROLLING_RESISTANCE_MOMENT);
        rolling_resistance_moment.clear();

        const auto& central_node = p_element->GetGeometry()[0];
        const double coeff_acc                            = central_node.FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) / dt;
        const array_1d<double, 3>& ang_velocity           = central_node.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        const array_1d<double, 3> initial_rotation_moment = coeff_acc * ang_velocity; // the moment needed to stop the spin in one time step

        const double MaxRotaMoment[3] = {initial_rotation_moment[0] + mContactMoment[0], initial_rotation_moment[1] + mContactMoment[1], initial_rotation_moment[2] + mContactMoment[2]};
        double CoordSystemMoment[3]   = {0.0};

        double max_rota_moment_modulus_inv = 1.0 / DEM_MODULUS_3(MaxRotaMoment);
        CoordSystemMoment[0]  = MaxRotaMoment[0] * max_rota_moment_modulus_inv;
        CoordSystemMoment[1]  = MaxRotaMoment[1] * max_rota_moment_modulus_inv;
        CoordSystemMoment[2]  = MaxRotaMoment[2] * max_rota_moment_modulus_inv;

        const double MR_now = DEM_INNER_PRODUCT_3(CoordSystemMoment, CoordSystemMoment) * mRollingResistance * mRollingResistance;
        const double MR_max = DEM_INNER_PRODUCT_3(MaxRotaMoment, MaxRotaMoment);

        if (MR_max > MR_now) {
            mContactMoment[0] -= CoordSystemMoment[0] * mRollingResistance;
            mContactMoment[1] -= CoordSystemMoment[1] * mRollingResistance;
            mContactMoment[2] -= CoordSystemMoment[2] * mRollingResistance;

            rolling_resistance_moment[0] -= CoordSystemMoment[0] * mRollingResistance;
            rolling_resistance_moment[1] -= CoordSystemMoment[1] * mRollingResistance;
            rolling_resistance_moment[2] -= CoordSystemMoment[2] * mRollingResistance;
        }
        else {
            rolling_resistance_moment = - mContactMoment;
            mContactMoment = - initial_rotation_moment;
        }

        KRATOS_CATCH("")
    }

}//namespace Kratos