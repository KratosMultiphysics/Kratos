// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019
// Saffman (1965, 1968)

#if !defined(SDEM_SAFFMAN_LIFT_LAW_H_INCLUDED)
#define SDEM_SAFFMAN_LIFT_LAW_H_INCLUDED

#include "vorticity_induced_lift_law.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) SaffmanLiftLaw : public VorticityInducedLiftLaw {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(SaffmanLiftLaw);

        // TODO: make mDoApplyFaxenCorrections an option
        SaffmanLiftLaw(){}

        SaffmanLiftLaw(Parameters r_parameters);

        ~SaffmanLiftLaw(){}

        VorticityInducedLiftLaw::Pointer Clone() const override;

        void Initialize(const ProcessInfo& r_process_info) override;

        std::string GetTypeOfLaw() override;

        void ComputeForce(Geometry<Node >& r_geometry,
                          const double reynolds_number,
                          double particle_radius,
                          double fluid_density,
                          double fluid_kinematic_viscosity,
                          array_1d<double, 3>& minus_slip_velocity,
                          array_1d<double, 3>& vorticity_induced_lift,
                          const ProcessInfo& r_current_process_info) override;

    private:

        friend class Serializer;
        double ComputeSaffmanLiftCoefficient(const double fluid_density,
                                             const double fluid_kinematic_viscosity,
                                             const double particle_radius,
                                             const double norm_of_vorticity);

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, VorticityInducedLiftLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, VorticityInducedLiftLaw)
        }

    }; //class SaffmanLiftLaw

} // Namespace Kratos

#endif /* SDEM_SAFFMAN_LIFT_LAW_H_INCLUDED  defined */
