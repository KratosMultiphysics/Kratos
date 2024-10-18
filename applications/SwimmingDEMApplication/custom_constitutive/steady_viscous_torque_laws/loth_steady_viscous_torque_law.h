// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019
// Loth (2008)
// Valid for Re_rot < 2000, Re_p < 20

#if !defined(SDEM_LOTH_STEADY_VISCOUS_TORQUE_LAW_H_INCLUDED)
#define SDEM_LOTH_STEADY_VISCOUS_TORQUE_LAW_H_INCLUDED

#include "rubinow_and_keller_torque_law.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) LothSteadyViscousTorqueLaw : public RubinowAndKellerTorqueLaw {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(LothSteadyViscousTorqueLaw);

        // TODO: make mDoApplyFaxenCorrections an option
        LothSteadyViscousTorqueLaw(){}

        LothSteadyViscousTorqueLaw(Parameters r_parameters);

        ~LothSteadyViscousTorqueLaw(){}

        SteadyViscousTorqueLaw::Pointer Clone() const override;

        void Initialize(const ProcessInfo& r_process_info) override;

        std::string GetTypeOfLaw() override;

        void ComputeMoment(Geometry<Node >& r_geometry,
                           const double reynolds_number,
                           double particle_radius,
                           double fluid_density,
                           double fluid_kinematic_viscosity,
                           array_1d<double, 3>& minus_slip_velocity,
                           array_1d<double, 3>& rotation_induced_lift,
                           const ProcessInfo& r_current_process_info) override;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SteadyViscousTorqueLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SteadyViscousTorqueLaw)
        }

    }; //class LothSteadyViscousTorqueLaw

} // Namespace Kratos

#endif /* SDEM_LOTH_STEADY_VISCOUS_TORQUE_LAW_H_INCLUDED  defined */
