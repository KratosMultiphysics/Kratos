// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#if !defined(SDEM_VIRTUAL_MASS_FORCE_LAW_H_INCLUDED)
#define SDEM_VIRTUAL_MASS_FORCE_LAW_H_INCLUDED

#include <string>
#include <iostream>
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/model_part.h"
#include "containers/flags.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) VirtualMassForceLaw : public Flags {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(VirtualMassForceLaw);

        VirtualMassForceLaw(): mLastVirtualMassAddedMass(0.0){}

        VirtualMassForceLaw(Parameters r_parameters): mLastVirtualMassAddedMass(0.0){}

        ~VirtualMassForceLaw(){}

        virtual VirtualMassForceLaw::Pointer Clone() const;

        virtual void Initialize(const ProcessInfo& r_process_info);

        void SetVirtualMassForceLawInProperties(Properties::Pointer pProp) const;

        virtual std::string GetTypeOfLaw();

        virtual void ComputeForce(Geometry<Node >& r_geometry,
                                  const double fluid_density,
                                  const double displaced_volume,
                                  array_1d<double, 3>& virtual_mass_flow_force,
                                  const ProcessInfo& r_current_process_info){}

        virtual double GetAddedMass(Geometry<Node >& r_geometry,
                                    double fluid_density,
                                    const ProcessInfo& r_current_process_info){return mLastVirtualMassAddedMass;}

        double ComputeParticleAccelerationNumber(const double particle_radius,
                                                 const array_1d<double, 3>& minus_slip_velocity,
                                                 const array_1d<double, 3>& minus_slip_acceleration);

    protected:

        double mLastVirtualMassAddedMass;

        virtual double GetVirtualMassCoefficient(Geometry<Node >& r_geometry,
                                                 const array_1d<double, 3>& minus_slip_acc){return 0.0;}

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
        }

    }; //class VirtualMassForceLaw

KRATOS_DEFINE_APPLICATION_VARIABLE(SWIMMING_DEM_APPLICATION, VirtualMassForceLaw::Pointer, SDEM_VIRTUAL_MASS_FORCE_LAW_POINTER)

} // Namespace Kratos

#endif /* SDEM_VIRTUAL_MASS_FORCE_LAW_H_INCLUDED  defined */
