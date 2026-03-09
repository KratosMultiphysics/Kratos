#include "power_law_hydrodynamic_interaction_law.h"
#include "swimming_DEM_application.h"

namespace Kratos {

    void PowerLawFluidHydrodynamicInteractionLaw::Initialize(const ProcessInfo& r_process_info)
    {}

    std::string PowerLawFluidHydrodynamicInteractionLaw::GetTypeOfLaw() {
        std::string type_of_law = "PowerLawFluidHydrodynamicInteractionLaw";
        return type_of_law;
    }

    PowerLawFluidHydrodynamicInteractionLaw::Pointer PowerLawFluidHydrodynamicInteractionLaw::Clone() const {
        PowerLawFluidHydrodynamicInteractionLaw::Pointer p_clone(new PowerLawFluidHydrodynamicInteractionLaw(*this));
        return p_clone;
    }

    PowerLawFluidHydrodynamicInteractionLaw::~PowerLawFluidHydrodynamicInteractionLaw(){}

    double PowerLawFluidHydrodynamicInteractionLaw::ComputeShahParticleReynoldsNumber(const double particle_radius,
                                                                                      const double fluid_density,
                                                                                      const double consistency_index,
                                                                                      const double flow_behavior_index,
                                                                                      const double modulus_of_minus_slip_velocity)
    {
        // This function is consistent with Shah 2007 (doi:10.1016/j.ijmultiphaseﬂow.2006.06.006)
        // int coefficient = use_max_shear_rate ? 3 : 2;
        const double& K = consistency_index;
        const double& n = flow_behavior_index;
        return 2 * std::pow(particle_radius, n) * std::pow(modulus_of_minus_slip_velocity, 2 - n) * fluid_density / K;
    }

    void PowerLawFluidHydrodynamicInteractionLaw::ComputeDragForce(SphericParticle* p_particle,
                                                                   double particle_radius,
                                                                   double fluid_density,
                                                                   double fluid_kinematic_viscosity,
                                                                   array_1d<double, 3>& minus_slip_velocity,
                                                                   array_1d<double, 3>& drag_force,
                                                                   const ProcessInfo& r_current_process_info)
    {
        const double consistency_index = r_current_process_info[POWER_LAW_K];
        const double flow_behavior_index = r_current_process_info[POWER_LAW_N];
        const double reynolds_number = ComputeShahParticleReynoldsNumber(particle_radius,
                                                                         fluid_density,
                                                                         consistency_index,
                                                                         flow_behavior_index,
                                                                         SWIMMING_MODULUS_3(minus_slip_velocity));
        mpDragLaw->ComputeForce(p_particle,
                                reynolds_number,
                                particle_radius,
                                fluid_density,
                                fluid_kinematic_viscosity,
                                minus_slip_velocity,
                                drag_force,
                                r_current_process_info);

    }

} // Namespace Kratos
