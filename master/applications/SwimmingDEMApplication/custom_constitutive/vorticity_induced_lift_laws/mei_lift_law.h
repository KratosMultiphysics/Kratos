// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019
// Mei (1992)

#if !defined(SDEM_MEI_LIFT_LAW_H_INCLUDED)
#define SDEM_MEI_LIFT_LAW_H_INCLUDED

#include "saffman_lift_law.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) MeiLiftLaw : public SaffmanLiftLaw {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(MeiLiftLaw);

        // TODO: make mDoApplyFaxenCorrections an option
        MeiLiftLaw(){}

        MeiLiftLaw(Parameters r_parameters);

        ~MeiLiftLaw(){}

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

        double ComputeShearReynoldsNumber(const double particle_radius,
                                          const double fluid_kinematic_viscosity,
                                          const double norm_of_vorticity);

        double ComputeMeiCorrectionOnSaffmanCoefficient(const double reynolds_number,
                                                        const double fluid_kinematic_viscosity,
                                                        const double particle_radius,
                                                        const double norm_of_vorticity);

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, VorticityInducedLiftLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, VorticityInducedLiftLaw)
        }

    }; //class MeiLiftLaw

} // Namespace Kratos

#endif /* SDEM_MEI_LIFT_LAW_H_INCLUDED  defined */
