// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#if !defined(SDEM_ZUBER_VIRTUAL_MASS_FORCE_LAW_H_INCLUDED)
#define SDEM_ZUBER_VIRTUAL_MASS_FORCE_LAW_H_INCLUDED

#include "auton_hunt_prudhomme_virtual_mass_force_law.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) ZuberVirtualMassForceLaw : public AutonHuntPrudhommeVirtualMassForceLaw {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(ZuberVirtualMassForceLaw);

        // TODO: make mDoApplyFaxenCorrections an option
        ZuberVirtualMassForceLaw(){}

        ZuberVirtualMassForceLaw(Parameters r_parameters);

        ~ZuberVirtualMassForceLaw(){}

        VirtualMassForceLaw::Pointer Clone() const override;

        void Initialize(const ProcessInfo& r_process_info) override;

        std::string GetTypeOfLaw() override;

    protected:

        double GetVirtualMassCoefficient(Geometry<Node >& r_geometry,
                                         const array_1d<double, 3>& minus_slip_acc) override;
    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, VirtualMassForceLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, VirtualMassForceLaw)
        }

    }; //class ZuberVirtualMassForceLaw

} // Namespace Kratos

#endif /* SDEM_ZUBER_VIRTUAL_MASS_FORCE_LAW_H_INCLUDED  defined */
