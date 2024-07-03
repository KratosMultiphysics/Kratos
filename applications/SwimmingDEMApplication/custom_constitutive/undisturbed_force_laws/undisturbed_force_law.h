// Author: Joaquin Gonzalez-Usua (jgonzalez@cimne.upc.edu)
// Date: January 2024

#if !defined(SDEM_UNDISTURBED_FORCE_LAW_H_INCLUDED)
#define SDEM_UNDISTURBED_FORCE_LAW_H_INCLUDED

#include <string>
#include <iostream>
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/model_part.h"
#include "containers/flags.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) UndisturbedForceLaw : public Flags {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(UndisturbedForceLaw);

        UndisturbedForceLaw(){}

        UndisturbedForceLaw(Parameters r_parameters);

        ~UndisturbedForceLaw(){}

        virtual UndisturbedForceLaw::Pointer Clone() const;

        virtual void Initialize(const ProcessInfo& r_process_info);

        void SetUndisturbedForceLawInProperties(Properties::Pointer pProp) const;

        virtual std::string GetTypeOfLaw();

        virtual void ComputeForce(Geometry<Node >& r_geometry,
                                  const double fluid_density,
                                  const double displaced_volume,
                                  array_1d<double, 3>& undisturbed_flow_force,
                                  const ProcessInfo& r_current_process_info);

    protected:

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
        }

    }; //class UndisturbedForceLaw


} // Namespace Kratos

#endif /* SDEM_UNDISTURBED_FORCE_LAW_H_INCLUDED  defined */