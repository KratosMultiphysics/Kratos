// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#if !defined(SDEM_VORTICITY_INDUCED_LIFT_LAW_H_INCLUDED)
#define SDEM_VORTICITY_INDUCED_LIFT_LAW_H_INCLUDED

#include <string>
#include <iostream>
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/model_part.h"
#include "containers/flags.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) VorticityInducedLiftLaw : public Flags {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(VorticityInducedLiftLaw);

        VorticityInducedLiftLaw(){}

        VorticityInducedLiftLaw(Parameters r_parameters){}

        ~VorticityInducedLiftLaw(){}

        virtual VorticityInducedLiftLaw::Pointer Clone() const;

        virtual void Initialize(const ProcessInfo& r_process_info);

        void SetVorticityInducedLiftLawInProperties(Properties::Pointer pProp) const;

        virtual std::string GetTypeOfLaw();

        virtual void ComputeForce(Geometry<Node >& r_geometry,
                                  const double reynolds_number,
                                  double particle_radius,
                                  double fluid_density,
                                  double fluid_kinematic_viscosity,
                                  array_1d<double, 3>& minus_slip_velocity,
                                  array_1d<double, 3>& vorticity_induced_lift,
                                  const ProcessInfo& r_current_process_info){}


    protected:

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
        }

    }; //class VorticityInducedLiftLaw

KRATOS_DEFINE_APPLICATION_VARIABLE(SWIMMING_DEM_APPLICATION, VorticityInducedLiftLaw::Pointer, SDEM_VORTICITY_INDUCED_LIFT_LAW_POINTER)

} // Namespace Kratos

#endif /* SDEM_VORTICITY_INDUCED_LIFT_LAW_H_INCLUDED  defined */
