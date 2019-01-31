#include "hydrodynamic_interaction_law.h"
#include "swimming_DEM_application.h"

namespace Kratos {

    HydrodynamicInteractionLaw::HydrodynamicInteractionLaw(const DragLaw& r_drag_law):
        mpDragLaw(r_drag_law.Clone()) {}

    HydrodynamicInteractionLaw::HydrodynamicInteractionLaw(const HydrodynamicInteractionLaw &rHydrodynamicInteractionLaw)
    {
        mpDragLaw = rHydrodynamicInteractionLaw.CloneDragLaw();
    }

    void HydrodynamicInteractionLaw::Initialize(const ProcessInfo& r_process_info) {
    }

    void HydrodynamicInteractionLaw::SetHydrodynamicInteractionLawInProperties(Properties::Pointer pProp, bool verbose) const {
        pProp->SetValue(SDEM_HYDRODYNAMIC_INTERACTION_LAW_POINTER, this->Clone());
    }

    std::string HydrodynamicInteractionLaw::GetTypeOfLaw() {
        std::string type_of_law = "HydrodynamicInteractionLaw";
        return type_of_law;
    }

    HydrodynamicInteractionLaw::Pointer HydrodynamicInteractionLaw::Clone() const {
        HydrodynamicInteractionLaw::Pointer p_clone(new HydrodynamicInteractionLaw(*this));
        return p_clone;
    }

    DragLaw::Pointer HydrodynamicInteractionLaw::CloneDragLaw() const {
        return mpDragLaw->Clone();
    }

    HydrodynamicInteractionLaw::~HydrodynamicInteractionLaw(){}

    double HydrodynamicInteractionLaw::ComputeParticleReynoldsNumber(const double particle_radius,
                                                                     const double fluid_kinematic_viscosity,
                                                                     const double modulus_of_slip_velocity)
    {
        return 2 * particle_radius * modulus_of_slip_velocity / fluid_kinematic_viscosity;
    }

    void HydrodynamicInteractionLaw::ComputeDragForce(Geometry<Node<3> >& r_geometry,
                                                      double particle_radius,
                                                      double fluid_density,
                                                      double fluid_kinematic_viscosity,
                                                      array_1d<double, 3>& slip_velocity,
                                                      array_1d<double, 3>& drag_force,
                                                      const ProcessInfo& r_current_process_info)
    {
        const double reynolds_number = ComputeParticleReynoldsNumber(particle_radius,
                                                                     fluid_kinematic_viscosity,
                                                                     SWIMMING_MODULUS_3(slip_velocity));
        mpDragLaw->ComputeForce(r_geometry,
                                reynolds_number,
                                particle_radius,
                                fluid_density,
                                fluid_kinematic_viscosity,
                                slip_velocity,
                                drag_force,
                                r_current_process_info);

    }

} // Namespace Kratos
