#ifndef KRATOS_SDEM_HYDRODYNAMIC_INTERACTION_LAW_H
#define KRATOS_SDEM_HYDRODYNAMIC_INTERACTION_LAW_H

#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/model_part.h"
#include "containers/flags.h"

#include "buoyancy_laws/buoyancy_law.h"
#include "drag_laws/drag_law.h"
#include "inviscid_force_laws/inviscid_force_law.h"
#include "history_force_laws/history_force_law.h"
#include "vorticity_induced_lift_laws/vorticity_induced_lift_law.h"
#include "rotation_induced_lift_laws/rotation_induced_lift_law.h"
#include "steady_viscous_torque_laws/steady_viscous_torque_law.h"
#include "../DEMApplication/custom_elements/spheric_particle.h"

namespace Kratos {

class KRATOS_API(SWIMMING_DEM_APPLICATION) HydrodynamicInteractionLaw : public Flags {

public:
    typedef Node <3> NodeType;

    // Pointer types for HydrodynamicInteractionLaw
    KRATOS_CLASS_POINTER_DEFINITION(HydrodynamicInteractionLaw);
    HydrodynamicInteractionLaw(){}
    HydrodynamicInteractionLaw(Properties::Pointer pProp, Parameters& r_hydrodynamic_parameters);
    HydrodynamicInteractionLaw(const HydrodynamicInteractionLaw &rHydrodynamicInteractionLaw);

    void SetBuoyancyLaw(const BuoyancyLaw& r_law){mpBuoyancyLaw = r_law.Clone();}
    void SetDragLaw(const DragLaw& r_law){mpDragLaw = r_law.Clone();}
    void SetInviscidForceLaw(const InviscidForceLaw& r_law){mpInviscidForceLaw = r_law.Clone();}
    void SetHistoryForceLaw(const HistoryForceLaw& r_law){mpHistoryForceLaw = r_law.Clone();}
    void SetVorticityInducedLiftLaw(const VorticityInducedLiftLaw& r_law){mpVorticityInducedLiftLaw = r_law.Clone();}
    void SetRotationInducedLiftLaw(const RotationInducedLiftLaw& r_law){mpRotationInducedLiftLaw = r_law.Clone();}
    void SetSteadyViscousTorqueLaw(const SteadyViscousTorqueLaw& r_law){mpSteadyViscousTorqueLaw = r_law.Clone();}

    virtual void Initialize(const ProcessInfo& r_process_info);

    virtual void SetHydrodynamicInteractionLawInProperties(Properties::Pointer pProp, bool verbose = true) const;

    virtual std::string GetTypeOfLaw();

    /// Destructor

    virtual ~HydrodynamicInteractionLaw();

    HydrodynamicInteractionLaw::Pointer Clone() const;

    virtual BuoyancyLaw::Pointer CloneBuoyancyLaw() const;
    virtual DragLaw::Pointer CloneDragLaw() const;
    virtual InviscidForceLaw::Pointer CloneInviscidForceLaw() const;
    virtual HistoryForceLaw::Pointer CloneHistoryForceLaw() const;
    virtual VorticityInducedLiftLaw::Pointer CloneVorticityInducedLiftLaw() const;
    virtual RotationInducedLiftLaw::Pointer CloneRotationInducedLiftLaw() const;
    virtual SteadyViscousTorqueLaw::Pointer CloneSteadyViscousTorqueLaw() const;

    double ComputeParticleReynoldsNumber(const double particle_radius,
                                         const double fluid_kinematic_viscosity,
                                         const double modulus_of_minus_slip_velocity);

    virtual void ComputeBuoyancyForce(Geometry<Node<3> >& r_geometry,
                                      const double fluid_density,
                                      const double displaced_volume,
                                      const array_1d<double, 3>& body_force,
                                      array_1d<double, 3>& buoyancy,
                                      const ProcessInfo& r_current_process_info);

    void ComputeDragForce(SphericParticle* p_particle,
                                  double particle_radius,
                                  double fluid_density,
                                  double fluid_kinematic_viscosity,
                                  array_1d<double, 3>& minus_slip_velocity,
                                  array_1d<double, 3>& drag_force,
                                  const ProcessInfo& r_current_process_info);

    virtual void ComputeInviscidForce(Geometry<Node<3> >& r_geometry,
                                      const double fluid_density,
                                      const double displaced_volume,
                                      array_1d<double, 3>& virtual_mass_plus_undisturbed_flow_force,
                                      const ProcessInfo& r_current_process_info);

    virtual double GetInviscidAddedMass(Geometry<Node<3> >& r_geometry,
                                        double fluid_density,
                                        const ProcessInfo& r_current_process_info);

    virtual void ComputeHistoryForce(Geometry<Node<3> >& r_geometry,
                                     double particle_radius,
                                     double fluid_density,
                                     double fluid_kinematic_viscosity,
                                     array_1d<double, 3>& minus_slip_velocity,
                                     array_1d<double, 3>& drag_force,
                                     const ProcessInfo& r_current_process_info);

    virtual double GetHistoryForceAddedMass(Geometry<Node<3> >& r_geometry,
                                            const ProcessInfo& r_current_process_info);

    virtual void ComputeVorticityInducedLift(Geometry<Node<3> >& r_geometry,
                                             double particle_radius,
                                             double fluid_density,
                                             double fluid_kinematic_viscosity,
                                             array_1d<double, 3>& minus_slip_velocity,
                                             array_1d<double, 3>& vorticity_induced_lift,
                                             const ProcessInfo& r_current_process_info);

    virtual void ComputeRotationInducedLift(Geometry<Node<3> >& r_geometry,
                                            double particle_radius,
                                            double fluid_density,
                                            double fluid_kinematic_viscosity,
                                            array_1d<double, 3>& minus_slip_velocity,
                                            array_1d<double, 3>& rotation_induced_lift,
                                            const ProcessInfo& r_current_process_info);

    virtual void ComputeSteadyViscousTorque(Geometry<Node<3> >& r_geometry,
                                            double particle_radius,
                                            double fluid_density,
                                            double fluid_kinematic_viscosity,
                                            array_1d<double, 3>& minus_slip_velocity,
                                            array_1d<double, 3>& steady_viscous_torque,
                                            const ProcessInfo& r_current_process_info);

protected:
    BuoyancyLaw::Pointer mpBuoyancyLaw;
    DragLaw::Pointer mpDragLaw;
    InviscidForceLaw::Pointer mpInviscidForceLaw;
    HistoryForceLaw::Pointer mpHistoryForceLaw;
    VorticityInducedLiftLaw::Pointer mpVorticityInducedLiftLaw;
    RotationInducedLiftLaw::Pointer mpRotationInducedLiftLaw;
    SteadyViscousTorqueLaw::Pointer mpSteadyViscousTorqueLaw;

private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)

    }

    virtual void load(Serializer& rSerializer) override {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
    }

}; // Class HydrodynamicInteractionLaw : public MainCL

KRATOS_DEFINE_APPLICATION_VARIABLE(SWIMMING_DEM_APPLICATION, HydrodynamicInteractionLaw::Pointer, SDEM_HYDRODYNAMIC_INTERACTION_LAW_POINTER)


} // Namespace Kratos

#endif // KRATOS_SDEM_HYDRODYNAMIC_INTERACTION_LAW_H
