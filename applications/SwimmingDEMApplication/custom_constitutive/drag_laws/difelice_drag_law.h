// Author: Joaquin Gonzalez-Usua (jgonzalez@cimne.upc.edu)
// Date: May 2023

#if !defined(SDEM_DIFELICE_DRAG_LAW_H_INCLUDED)
#define SDEM_DIFELICE_DRAG_LAW_H_INCLUDED
#include <algorithm>
#include "stokes_drag_law.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) DiFeliceDragLaw : public StokesDragLaw {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(DiFeliceDragLaw);

        DiFeliceDragLaw(): StokesDragLaw(){}

        DiFeliceDragLaw(Parameters& r_parameters): StokesDragLaw(r_parameters){}

        ~DiFeliceDragLaw(){}

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

        array_1d<double, 3> CalculateParticleVelocityMean(SphericParticle* p_particle);

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DragLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DragLaw)
        }

    }; //class DiFeliceDragLaw

} // Namespace Kratos

#endif /* SDEM_DIFELICE_DRAG_LAW_H_INCLUDED  defined */