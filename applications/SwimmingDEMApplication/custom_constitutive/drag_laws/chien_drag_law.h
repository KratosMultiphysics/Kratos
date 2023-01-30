#if !defined(SDEM_CHIEN_DRAG_LAW_H_INCLUDED)
#define SDEM_CHIEN_DRAG_LAW_H_INCLUDED

#include "stokes_drag_law.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) ChienDragLaw : public StokesDragLaw {

    public:
        typedef Node <3> NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(ChienDragLaw);

        ChienDragLaw(): StokesDragLaw(){}

        ChienDragLaw(Parameters r_parameters): StokesDragLaw(r_parameters){}

        ~ChienDragLaw(){}

        DragLaw::Pointer Clone() const override;

        void Initialize(const ProcessInfo& r_process_info) override;

        std::string GetTypeOfLaw() override;

        void ComputeForce(SphericParticle* r_geometry,
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

    }; //class ChienDragLaw

} // Namespace Kratos

#endif /* SDEM_CHIEN_DRAG_LAW_H_INCLUDED  defined */
