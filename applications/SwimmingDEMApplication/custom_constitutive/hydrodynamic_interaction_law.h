#ifndef KRATOS_SDEM_HYDRODYNAMIC_INTERACTION_LAW_H
#define KRATOS_SDEM_HYDRODYNAMIC_INTERACTION_LAW_H

#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/model_part.h"
#include "containers/flags.h"

#include "drag_law.h"

namespace Kratos {

class KRATOS_API(SWIMMING_DEM_APPLICATION) HydrodynamicInteractionLaw : public Flags {

public:
    typedef Node <3> NodeType;

    // Pointer types for HydrodynamicInteractionLaw
    KRATOS_CLASS_POINTER_DEFINITION(HydrodynamicInteractionLaw);

    HydrodynamicInteractionLaw();

    HydrodynamicInteractionLaw(const HydrodynamicInteractionLaw &rHydrodynamicInteractionLaw);

    virtual void Initialize(const ProcessInfo& r_process_info);

    virtual void SetHydrodynamicInteractionLawInProperties(Properties::Pointer pProp, bool verbose = true) const;

    virtual std::string GetTypeOfLaw();

    /// Destructor

    virtual ~HydrodynamicInteractionLaw();

    virtual HydrodynamicInteractionLaw::Pointer Clone() const;

    void ComputeDragForce(NodeType& node,
                            double particle_radius,
                            double fluid_density,
                            double fluid_kinematic_viscosity,
                            array_1d<double, 3>& slip_velocity,
                            array_1d<double, 3>& drag_force,
                            const ProcessInfo& r_current_process_info);

protected:
    DragLaw::Pointer mpDragLaw;

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
