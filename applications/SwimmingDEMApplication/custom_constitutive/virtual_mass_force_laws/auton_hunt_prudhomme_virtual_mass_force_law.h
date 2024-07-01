// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#if !defined(SDEM_STANDARD_VIRTUAL_MASS_FORCE_LAW_H_INCLUDED)
#define SDEM_STANDARD_VIRTUAL_MASS_FORCE_LAW_H_INCLUDED

#include "virtual_mass_force_law.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) AutonHuntPrudhommeVirtualMassForceLaw : public VirtualMassForceLaw {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(AutonHuntPrudhommeVirtualMassForceLaw);

        // TODO: make mDoApplyFaxenCorrections an option
        AutonHuntPrudhommeVirtualMassForceLaw(){}

        AutonHuntPrudhommeVirtualMassForceLaw(Parameters r_parameters);

        ~AutonHuntPrudhommeVirtualMassForceLaw(){}

        VirtualMassForceLaw::Pointer Clone() const override;

        void Initialize(const ProcessInfo& r_process_info) override;

        std::string GetTypeOfLaw() override;

        void ComputeForce(Geometry<Node >& r_geometry,
                          const double fluid_density,
                          const double displaced_volume,
                          array_1d<double, 3>& virtual_mass_force,
                          const ProcessInfo& r_current_process_info) override;

    protected:

        double GetVirtualMassCoefficient(Geometry<Node >& r_geometry,
                                         const array_1d<double, 3>& minus_slip_acc) override;

    private:
        bool mDoApplyFaxenCorrections;

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, VirtualMassForceLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, VirtualMassForceLaw)
        }

    }; //class AutonHuntPrudhommeVirtualMassForceLaw

} // Namespace Kratos

#endif /* SDEM_STANDARD_VIRTUAL_MASS_FORCE_LAW_H_INCLUDED  defined */
