#if !defined(DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_2D_H_INCLUDED)
#define DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_2D_H_INCLUDED

/* Project includes */
#include "DEM_KDEM_with_damage_parallel_bond_CL.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_KDEM_with_damage_parallel_bond_2D : public DEM_KDEM_with_damage_parallel_bond {

        //typedef DEM_KDEM_with_damage_parallel_bond BaseClassType;

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM_with_damage_parallel_bond_2D);

        DEM_KDEM_with_damage_parallel_bond_2D() {}

        ~DEM_KDEM_with_damage_parallel_bond_2D() {}

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        void CalculateContactArea(double radius, double other_radius, double& calculation_area) override;


    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

} // namespace Kratos
#endif /* DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_H_INCLUDED defined */
