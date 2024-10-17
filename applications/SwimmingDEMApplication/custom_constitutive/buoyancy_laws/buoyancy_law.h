// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#if !defined(SDEM_BUOYANCY_FORCE_FORCE_LAW_H_INCLUDED)
#define SDEM_BUOYANCY_FORCE_FORCE_LAW_H_INCLUDED

#include <string>
#include <iostream>
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/model_part.h"
#include "containers/flags.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) BuoyancyLaw : public Flags {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(BuoyancyLaw);

        BuoyancyLaw(){}

        BuoyancyLaw(Parameters r_hydrodynamic_parameters){}

        ~BuoyancyLaw(){}

        virtual BuoyancyLaw::Pointer Clone() const;

        virtual void Initialize(const ProcessInfo& r_process_info);

        void SetBuoyancyLawInProperties(Properties::Pointer pProp) const;

        virtual std::string GetTypeOfLaw();

        virtual void ComputeForce(Geometry<Node >& r_geometry,
                                  const double fluid_density,
                                  const double displaced_volume,
                                  const array_1d<double, 3>& body_force,
                                  array_1d<double, 3>& buoyancy,
                                  const ProcessInfo& r_current_process_info) {}

    protected:

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
        }

    }; //class BuoyancyLaw

KRATOS_DEFINE_APPLICATION_VARIABLE(SWIMMING_DEM_APPLICATION, BuoyancyLaw::Pointer, SDEM_BUOYANCY_LAW_POINTER)

} // Namespace Kratos

#endif /* SDEM_BUOYANCY_FORCE_FORCE_LAW_H_INCLUDED  defined */
