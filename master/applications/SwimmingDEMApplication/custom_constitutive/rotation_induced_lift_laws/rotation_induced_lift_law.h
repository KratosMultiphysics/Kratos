// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#if !defined(SDEM_ROTATION_INDUCED_LIFT_LAW_H_INCLUDED)
#define SDEM_ROTATION_INDUCED_LIFT_LAW_H_INCLUDED

#include <string>
#include <iostream>
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/model_part.h"
#include "containers/flags.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) RotationInducedLiftLaw : public Flags {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(RotationInducedLiftLaw);

        RotationInducedLiftLaw(){}

        RotationInducedLiftLaw(Parameters r_parameters){}

        ~RotationInducedLiftLaw(){}

        virtual RotationInducedLiftLaw::Pointer Clone() const;

        virtual void Initialize(const ProcessInfo& r_process_info);

        void SetRotationInducedLiftLawInProperties(Properties::Pointer pProp) const;

        virtual std::string GetTypeOfLaw();

        virtual void ComputeForce(Geometry<Node >& r_geometry,
                                  const double reynolds_number,
                                  double particle_radius,
                                  double fluid_density,
                                  double fluid_kinematic_viscosity,
                                  array_1d<double, 3>& minus_slip_velocity,
                                  array_1d<double, 3>& rotation_induced_lift,
                                  const ProcessInfo& r_current_process_info){}


    protected:
        double ComputeParticleRotationReynoldsNumber(const double norm_of_slip_rot,
                                                     const double particle_radius,
                                                     const double fluid_kinematic_viscosity);

        double ComputeNondimensionalRotVelocity(const double norm_of_slip_vel,
                                                const double norm_of_slip_rot,
                                                const double particle_radius,
                                                const double fluid_kinematic_viscosity);
    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
        }

    }; //class RotationInducedLiftLaw

KRATOS_DEFINE_APPLICATION_VARIABLE(SWIMMING_DEM_APPLICATION, RotationInducedLiftLaw::Pointer, SDEM_ROTATION_INDUCED_LIFT_LAW_POINTER)

} // Namespace Kratos

#endif /* SDEM_ROTATION_INDUCED_LIFT_LAW_H_INCLUDED  defined */
