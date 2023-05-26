// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#if !defined(SDEM_DRAG_LAW_H_INCLUDED)
#define SDEM_DRAG_LAW_H_INCLUDED

#include <string>
#include <iostream>
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/model_part.h"
#include "containers/flags.h"
#include "../DEMApplication/custom_elements/spheric_particle.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) DragLaw : public Flags {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(DragLaw);

        DragLaw(){}

        DragLaw(Parameters r_parameters){}

        ~DragLaw(){}

        virtual DragLaw::Pointer Clone() const;

        virtual void Initialize(const ProcessInfo& r_process_info);

        void SetDragLawInProperties(Properties::Pointer pProp) const;

        virtual std::string GetTypeOfLaw();

        virtual void ComputeForce(SphericParticle* p_particle,
                                  const double reynolds_number,
                                  double particle_radius,
                                  double fluid_density,
                                  double fluid_kinematic_viscosity,
                                  array_1d<double, 3>& minus_slip_velocity,
                                  array_1d<double, 3>& drag_force,
                                  const ProcessInfo& r_current_process_info);

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
        }

    }; //class DragLaw

KRATOS_DEFINE_APPLICATION_VARIABLE(SWIMMING_DEM_APPLICATION, DragLaw::Pointer, SDEM_DRAG_LAW_POINTER)

} // Namespace Kratos

#endif /* SDEM_DRAG_LAW_H_INCLUDED  defined */
