// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#if !defined(SDEM_STANDARD_INVISCID_FORCE_LAW_H_INCLUDED)
#define SDEM_STANDARD_INVISCID_FORCE_LAW_H_INCLUDED

#include "inviscid_force_law.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) AutonHuntPrudhommeInviscidForceLaw : public InviscidForceLaw {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(AutonHuntPrudhommeInviscidForceLaw);

        // TODO: make mDoApplyFaxenCorrections an option
        AutonHuntPrudhommeInviscidForceLaw(){}

        AutonHuntPrudhommeInviscidForceLaw(Parameters r_parameters);

        ~AutonHuntPrudhommeInviscidForceLaw(){}

        InviscidForceLaw::Pointer Clone() const override;

        void Initialize(const ProcessInfo& r_process_info) override;

        std::string GetTypeOfLaw() override;

        void ComputeForce(Geometry<Node >& r_geometry,
                          const double fluid_density,
                          const double displaced_volume,
                          array_1d<double, 3>& virtual_mass_plus_undisturbed_flow_force,
                          const ProcessInfo& r_current_process_info) override;

    protected:

        double GetVirtualMassCoefficient(Geometry<Node >& r_geometry,
                                         const array_1d<double, 3>& minus_slip_acc) override;

    private:
        bool mDoApplyFaxenCorrections;

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, InviscidForceLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, InviscidForceLaw)
        }

    }; //class AutonHuntPrudhommeInviscidForceLaw

} // Namespace Kratos

#endif /* SDEM_STANDARD_INVISCID_FORCE_LAW_H_INCLUDED  defined */
