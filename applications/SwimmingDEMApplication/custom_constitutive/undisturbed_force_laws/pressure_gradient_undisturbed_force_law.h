// Author: Joaquin Gonzalez-Usua (jgonzalez@cimne.upc.edu)
// Date: January 2024

#if !defined(SDEM_PRESSURE_GRADIENT_UNDISTURBED_FORCE_LAW_H_INCLUDED)
#define SDEM_PRESSURE_GRADIENT_UNDISTURBED_FORCE_LAW_H_INCLUDED

#include "undisturbed_force_law.h"
//#include "auton_hunt_prudhomme_inviscid_force_law.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) PressureGradientUndisturbedForceLaw : public UndisturbedForceLaw {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(PressureGradientUndisturbedForceLaw);

        // TODO: make mDoApplyFaxenCorrections an option
        PressureGradientUndisturbedForceLaw(){}

        PressureGradientUndisturbedForceLaw(Parameters r_parameters);

        ~PressureGradientUndisturbedForceLaw(){}

        UndisturbedForceLaw::Pointer Clone() const override;

        void Initialize(const ProcessInfo& r_process_info) override;

        std::string GetTypeOfLaw() override;

        void ComputeForce(Geometry<Node >& r_geometry,
                        const double fluid_density,
                        const double displaced_volume,
                        array_1d<double, 3>& undisturbed_flow_force,
                        const ProcessInfo& r_current_process_info) override;

    protected:

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, UndisturbedForceLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, UndisturbedForceLaw)
        }

    }; //class PressureGradientUndisturbedForceLaw

} // Namespace Kratos

#endif /* SDEM_PRESSURE_GRADIENT_UNDISTURBED_FORCE_LAW_H_INCLUDED  defined */
