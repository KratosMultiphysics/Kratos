// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#if !defined(SDEM_ZUBER_INVISCID_FORCE_LAW_H_INCLUDED)
#define SDEM_ZUBER_INVISCID_FORCE_LAW_H_INCLUDED

#include "auton_hunt_prudhomme_inviscid_force_law.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) ZuberInviscidForceLaw : public AutonHuntPrudhommeInviscidForceLaw {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(ZuberInviscidForceLaw);

        // TODO: make mDoApplyFaxenCorrections an option
        ZuberInviscidForceLaw(){}

        ZuberInviscidForceLaw(Parameters r_parameters);

        ~ZuberInviscidForceLaw(){}

        InviscidForceLaw::Pointer Clone() const override;

        void Initialize(const ProcessInfo& r_process_info) override;

        std::string GetTypeOfLaw() override;

    protected:

        double GetVirtualMassCoefficient(Geometry<Node >& r_geometry,
                                         const array_1d<double, 3>& minus_slip_acc) override;
    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, InviscidForceLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, InviscidForceLaw)
        }

    }; //class ZuberInviscidForceLaw

} // Namespace Kratos

#endif /* SDEM_ZUBER_INVISCID_FORCE_LAW_H_INCLUDED  defined */
