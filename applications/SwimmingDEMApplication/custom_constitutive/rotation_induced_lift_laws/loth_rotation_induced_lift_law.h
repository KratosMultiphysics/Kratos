// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019
// Loth (2008)
// Valid for Re_p < 2000; nondimensional_minus_slip_rot_vel < 20

#if !defined(SDEM_LOTH_ROTATION_INDUCED_LIFT_LAW_H_INCLUDED)
#define SDEM_LOTH_ROTATION_INDUCED_LIFT_LAW_H_INCLUDED

#include "rubinow_and_keller_lift_law.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) LothRotationInducedLiftLaw : public RubinowAndKellerLiftLaw {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(LothRotationInducedLiftLaw);

        // TODO: make mDoApplyFaxenCorrections an option
        LothRotationInducedLiftLaw(){}

        LothRotationInducedLiftLaw(Parameters r_parameters);

        ~LothRotationInducedLiftLaw(){}

        RotationInducedLiftLaw::Pointer Clone() const override;

        void Initialize(const ProcessInfo& r_process_info) override;

        std::string GetTypeOfLaw() override;

        void ComputeForce(Geometry<Node >& r_geometry,
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
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, RotationInducedLiftLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, RotationInducedLiftLaw)
        }

    }; //class LothRotationInducedLiftLaw

} // Namespace Kratos

#endif /* SDEM_LOTH_ROTATION_INDUCED_LIFT_LAW_H_INCLUDED  defined */
