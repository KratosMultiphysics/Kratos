#include "swimming_DEM_application.h"
#include "hydrodynamic_interaction_law.h"

namespace Kratos {

    HydrodynamicInteractionLaw::HydrodynamicInteractionLaw(Properties::Pointer pProp, Parameters& r_hydrodynamic_parameters):
        mpBuoyancyLaw(BuoyancyLaw().Clone()),
        mpDragLaw(DragLaw().Clone()),
        mpInviscidForceLaw(InviscidForceLaw().Clone()),
        mpHistoryForceLaw(HistoryForceLaw().Clone()),
        mpVorticityInducedLiftLaw(VorticityInducedLiftLaw().Clone()),
        mpRotationInducedLiftLaw(RotationInducedLiftLaw().Clone()),
        mpSteadyViscousTorqueLaw(SteadyViscousTorqueLaw().Clone()){}

    HydrodynamicInteractionLaw::HydrodynamicInteractionLaw(const HydrodynamicInteractionLaw &rHydrodynamicInteractionLaw)
    {
        mpBuoyancyLaw = rHydrodynamicInteractionLaw.CloneBuoyancyLaw();
        mpDragLaw = rHydrodynamicInteractionLaw.CloneDragLaw();
        mpInviscidForceLaw = rHydrodynamicInteractionLaw.CloneInviscidForceLaw();
        mpHistoryForceLaw = rHydrodynamicInteractionLaw.CloneHistoryForceLaw();
        mpVorticityInducedLiftLaw = rHydrodynamicInteractionLaw.CloneVorticityInducedLiftLaw();
        mpRotationInducedLiftLaw = rHydrodynamicInteractionLaw.CloneRotationInducedLiftLaw();
        mpSteadyViscousTorqueLaw = rHydrodynamicInteractionLaw.CloneSteadyViscousTorqueLaw();
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

    BuoyancyLaw::Pointer HydrodynamicInteractionLaw::CloneBuoyancyLaw() const {
        return mpBuoyancyLaw->Clone();
    }

    DragLaw::Pointer HydrodynamicInteractionLaw::CloneDragLaw() const {
        return mpDragLaw->Clone();
    }

    InviscidForceLaw::Pointer HydrodynamicInteractionLaw::CloneInviscidForceLaw() const {
        return mpInviscidForceLaw->Clone();
    }

    HistoryForceLaw::Pointer HydrodynamicInteractionLaw::CloneHistoryForceLaw() const {
        return mpHistoryForceLaw->Clone();
    }

    VorticityInducedLiftLaw::Pointer HydrodynamicInteractionLaw::CloneVorticityInducedLiftLaw() const {
        return mpVorticityInducedLiftLaw->Clone();
    }

    RotationInducedLiftLaw::Pointer HydrodynamicInteractionLaw::CloneRotationInducedLiftLaw() const {
        return mpRotationInducedLiftLaw->Clone();
    }

    SteadyViscousTorqueLaw::Pointer HydrodynamicInteractionLaw::CloneSteadyViscousTorqueLaw() const {
        return mpSteadyViscousTorqueLaw->Clone();
    }

    HydrodynamicInteractionLaw::~HydrodynamicInteractionLaw(){}

    double HydrodynamicInteractionLaw::ComputeParticleReynoldsNumber(const double particle_radius,
                                                                     const double fluid_kinematic_viscosity,
                                                                     const double modulus_of_minus_slip_velocity)
    {
        return 2 * particle_radius * modulus_of_minus_slip_velocity / fluid_kinematic_viscosity;
    }

    void HydrodynamicInteractionLaw::ComputeBuoyancyForce(Geometry<Node<3> >& r_geometry,
                                                          const double fluid_density,
                                                          const double displaced_volume,
                                                          const array_1d<double, 3>& body_force,
                                                          array_1d<double, 3>& buoyancy,
                                                          const ProcessInfo& r_current_process_info)
    {
        mpBuoyancyLaw->ComputeForce(r_geometry,
                                    fluid_density,
                                    displaced_volume,
                                    body_force,
                                    buoyancy,
                                    r_current_process_info);
    }

    void HydrodynamicInteractionLaw::ComputeDragForce(SphericParticle* p_particle,
                                                      double particle_radius,
                                                      double fluid_density,
                                                      double fluid_kinematic_viscosity,
                                                      array_1d<double, 3>& minus_slip_velocity,
                                                      array_1d<double, 3>& drag_force,
                                                      const ProcessInfo& r_current_process_info)
    {
        const double reynolds_number = ComputeParticleReynoldsNumber(particle_radius,
                                                                     fluid_kinematic_viscosity,
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

    void HydrodynamicInteractionLaw::ComputeInviscidForce(Geometry<Node<3> >& r_geometry,
                                                          const double fluid_density,
                                                          const double displaced_volume,
                                                          array_1d<double, 3>& virtual_mass_plus_undisturbed_flow_force,
                                                          const ProcessInfo& r_current_process_info)
    {
        mpInviscidForceLaw->ComputeForce(r_geometry,
                                         fluid_density,
                                         displaced_volume,
                                         virtual_mass_plus_undisturbed_flow_force,
                                         r_current_process_info);
    }

    double HydrodynamicInteractionLaw::GetInviscidAddedMass(Geometry<Node<3> >& r_geometry,
                                                            double fluid_density,
                                                            const ProcessInfo& r_current_process_info)
    {
        return mpInviscidForceLaw->GetAddedMass(r_geometry,
                                                fluid_density,
                                                r_current_process_info);
    }

    void HydrodynamicInteractionLaw::ComputeHistoryForce(Geometry<Node<3> >& r_geometry,
                                                         double particle_radius,
                                                         double fluid_density,
                                                         double fluid_kinematic_viscosity,
                                                         array_1d<double, 3>& minus_slip_velocity,
                                                         array_1d<double, 3>& basset_force,
                                                         const ProcessInfo& r_current_process_info)
    {
        const double reynolds_number = ComputeParticleReynoldsNumber(particle_radius,
                                                                     fluid_kinematic_viscosity,
                                                                     SWIMMING_MODULUS_3(minus_slip_velocity));
        return mpHistoryForceLaw->ComputeForce(r_geometry,
                                               reynolds_number,
                                               particle_radius,
                                               fluid_density,
                                               fluid_kinematic_viscosity,
                                               minus_slip_velocity,
                                               basset_force,
                                               r_current_process_info);
    }

    double HydrodynamicInteractionLaw::GetHistoryForceAddedMass(Geometry<Node<3> >& r_geometry,
                                                                const ProcessInfo& r_current_process_info)
    {
        return mpHistoryForceLaw->GetAddedMass(r_geometry, r_current_process_info);
    }

    void HydrodynamicInteractionLaw::ComputeVorticityInducedLift(Geometry<Node<3> >& r_geometry,
                                                                 double particle_radius,
                                                                 double fluid_density,
                                                                 double fluid_kinematic_viscosity,
                                                                 array_1d<double, 3>& minus_slip_velocity,
                                                                 array_1d<double, 3>& vorticity_induced_lift,
                                                                 const ProcessInfo& r_current_process_info)
    {
        const double reynolds_number = ComputeParticleReynoldsNumber(particle_radius,
                                                                     fluid_kinematic_viscosity,
                                                                     SWIMMING_MODULUS_3(minus_slip_velocity));
        return mpVorticityInducedLiftLaw->ComputeForce(r_geometry,
                                                       reynolds_number,
                                                       particle_radius,
                                                       fluid_density,
                                                       fluid_kinematic_viscosity,
                                                       minus_slip_velocity,
                                                       vorticity_induced_lift,
                                                       r_current_process_info);
    }

    void HydrodynamicInteractionLaw::ComputeRotationInducedLift(Geometry<Node<3> >& r_geometry,
                                                                double particle_radius,
                                                                double fluid_density,
                                                                double fluid_kinematic_viscosity,
                                                                array_1d<double, 3>& minus_slip_velocity,
                                                                array_1d<double, 3>& rotation_induced_lift,
                                                                const ProcessInfo& r_current_process_info)
    {
        const double reynolds_number = ComputeParticleReynoldsNumber(particle_radius,
                                                                     fluid_kinematic_viscosity,
                                                                     SWIMMING_MODULUS_3(minus_slip_velocity));
        return mpRotationInducedLiftLaw->ComputeForce(r_geometry,
                                                      reynolds_number,
                                                      particle_radius,
                                                      fluid_density,
                                                      fluid_kinematic_viscosity,
                                                      minus_slip_velocity,
                                                      rotation_induced_lift,
                                                      r_current_process_info);
    }

    void HydrodynamicInteractionLaw::ComputeSteadyViscousTorque(Geometry<Node<3> >& r_geometry,
                                                                double particle_radius,
                                                                double fluid_density,
                                                                double fluid_kinematic_viscosity,
                                                                array_1d<double, 3>& minus_slip_velocity,
                                                                array_1d<double, 3>& steady_viscous_torque,
                                                                const ProcessInfo& r_current_process_info)
    {
        const double reynolds_number = ComputeParticleReynoldsNumber(particle_radius,
                                                                     fluid_kinematic_viscosity,
                                                                     SWIMMING_MODULUS_3(minus_slip_velocity));
        return mpSteadyViscousTorqueLaw->ComputeMoment(r_geometry,
                                                      reynolds_number,
                                                      particle_radius,
                                                      fluid_density,
                                                      fluid_kinematic_viscosity,
                                                      minus_slip_velocity,
                                                      steady_viscous_torque,
                                                      r_current_process_info);
    }

} // Namespace Kratos
