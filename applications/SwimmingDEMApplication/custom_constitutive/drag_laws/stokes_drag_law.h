// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#if !defined(SDEM_STOKES_DRAG_LAW_H_INCLUDED)
#define SDEM_STOKES_DRAG_LAW_H_INCLUDED

#include "drag_law.h"
#include "../DEMApplication/custom_elements/spheric_particle.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) StokesDragLaw : public DragLaw {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(StokesDragLaw);

        StokesDragLaw(){}

        StokesDragLaw(Parameters r_parameters){}

        ~StokesDragLaw(){}

        DragLaw::Pointer Clone() const override;

        void Initialize(const ProcessInfo& r_process_info) override;

        std::string GetTypeOfLaw() override;

        void ComputeForce(SphericParticle* p_particle,
                          const double reynolds_number,
                          double particle_radius,
                          double fluid_density,
                          double fluid_kinematic_viscosity,
                          array_1d<double, 3>& minus_slip_velocity,
                          array_1d<double, 3>& drag_force,
                          const ProcessInfo& r_current_process_info) override;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DragLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DragLaw)
        }

    }; //class StokesDragLaw

} // Namespace Kratos

#endif /* SDEM_STOKES_DRAG_LAW_H_INCLUDED  defined */
