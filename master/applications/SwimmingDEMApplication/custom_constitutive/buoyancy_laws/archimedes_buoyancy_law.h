// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#if !defined(SDEM_STANDARD_BUOYANCY_FORCE_LAW_H_INCLUDED)
#define SDEM_STANDARD_BUOYANCY_FORCE_LAW_H_INCLUDED

#include "buoyancy_law.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) ArchimedesBuoyancyLaw : public BuoyancyLaw {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(ArchimedesBuoyancyLaw);

        // TODO: make mDoApplyFaxenCorrections an option
        ArchimedesBuoyancyLaw(){}

        ArchimedesBuoyancyLaw(Parameters r_hydrodynamic_parameters){}

        ~ArchimedesBuoyancyLaw(){}

        BuoyancyLaw::Pointer Clone() const override;

        void Initialize(const ProcessInfo& r_process_info) override;

        std::string GetTypeOfLaw() override;

        void ComputeForce(Geometry<Node >& r_geometry,
                          const double fluid_density,
                          const double displaced_volume,
                          const array_1d<double, 3>& body_force,
                          array_1d<double, 3>& buoyancy,
                          const ProcessInfo& r_current_process_info) override;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BuoyancyLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BuoyancyLaw)
        }

    }; //class ArchimedesBuoyancyLaw

} // Namespace Kratos

#endif /* SDEM_STANDARD_BUOYANCY_FORCE_LAW_H_INCLUDED  defined */
