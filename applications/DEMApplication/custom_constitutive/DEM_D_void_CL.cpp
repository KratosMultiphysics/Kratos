/////////////////////////////////////////////////
// Author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu, chengshun.shang1996@gmail.com
// Date: June 2023
/////////////////////////////////////////////////

#include "DEM_D_void_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_void::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_void(*this));
        return p_clone;
    }

    std::unique_ptr<DEMDiscontinuumConstitutiveLaw> DEM_D_void::CloneUnique() {
        return Kratos::make_unique<DEM_D_void>();
    }

    std::string DEM_D_void::GetTypeOfLaw() {
        std::string type_of_law = "Void"; 
        return type_of_law;
    }

    void DEM_D_void::Check(Properties::Pointer pProp) const {}

    /////////////////////////
    // DEM-DEM INTERACTION //
    /////////////////////////

    void DEM_D_void::CalculateForces(const ProcessInfo& r_process_info,
                                                       const double OldLocalElasticContactForce[3],
                                                       double LocalElasticContactForce[3],
                                                       double LocalDeltDisp[3],
                                                       double LocalRelVel[3],
                                                       double indentation,
                                                       double previous_indentation,
                                                       double ViscoDampingLocalContactForce[3],
                                                       double& cohesive_force,
                                                       SphericParticle* element1,
                                                       SphericParticle* element2,
                                                       bool& sliding, double LocalCoordSystem[3][3]) {
        KRATOS_TRY
        LocalElasticContactForce[0]  = 0.0;
        LocalElasticContactForce[1]  = 0.0;
        LocalElasticContactForce[2]  = 0.0;
        ViscoDampingLocalContactForce[0] = 0.0;
        ViscoDampingLocalContactForce[1] = 0.0;
        ViscoDampingLocalContactForce[2] = 0.0;
        KRATOS_CATCH("")
    }

    /////////////////////////
    // DEM-FEM INTERACTION //
    /////////////////////////

    void DEM_D_void::CalculateForcesWithFEM(const ProcessInfo& r_process_info,
                                                            const double OldLocalElasticContactForce[3],
                                                            double LocalElasticContactForce[3],
                                                            double LocalDeltDisp[3],
                                                            double LocalRelVel[3],
                                                            double indentation,
                                                            double previous_indentation,
                                                            double ViscoDampingLocalContactForce[3],
                                                            double& cohesive_force,
                                                            SphericParticle* const element,
                                                            Condition* const wall,
                                                            bool& sliding) {
        KRATOS_TRY
        LocalElasticContactForce[0]  = 0.0;
        LocalElasticContactForce[1]  = 0.0;
        LocalElasticContactForce[2]  = 0.0;
        ViscoDampingLocalContactForce[0] = 0.0;
        ViscoDampingLocalContactForce[1] = 0.0;
        ViscoDampingLocalContactForce[2] = 0.0;
        KRATOS_CATCH("")
    }

} // namespace Kratos
