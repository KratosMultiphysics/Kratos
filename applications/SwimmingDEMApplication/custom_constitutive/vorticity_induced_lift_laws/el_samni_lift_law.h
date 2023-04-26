// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019
// El Samni, E.A. (1949), see paper by R. K. Clark (1994)

#if !defined(SDEM_EL_SAMNI_LIFT_LAW_H_INCLUDED)
#define SDEM_EL_SAMNI_LIFT_LAW_H_INCLUDED

#include "vorticity_induced_lift_law.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) ElSamniLiftLaw : public VorticityInducedLiftLaw {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(ElSamniLiftLaw);

        // TODO: make mDoApplyFaxenCorrections an option
        ElSamniLiftLaw(){}

        ElSamniLiftLaw(Parameters r_parameters);

        ~ElSamniLiftLaw(){}

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
        double ComputeElSamniLiftCoefficient(const double particle_radius,
                                             const double fluid_density,
                                             const double norm_of_slip_vel,
                                             const double vorticity_norm,
                                             const ProcessInfo& r_current_process_info);

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, VorticityInducedLiftLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, VorticityInducedLiftLaw)
        }

    }; //class ElSamniLiftLaw

} // Namespace Kratos

#endif /* SDEM_EL_SAMNI_LIFT_LAW_H_INCLUDED  defined */
