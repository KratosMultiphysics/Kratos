#if !defined(DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_CAPPED_H_INCLUDED)
#define DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_CAPPED_H_INCLUDED

/* Project includes */
#include "DEM_KDEM_with_damage_parallel_bond_CL.h"
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_KDEM_with_damage_parallel_bond_capped : public DEM_KDEM_with_damage_parallel_bond {

        typedef DEM_KDEM_with_damage_parallel_bond BaseClassType;

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM_with_damage_parallel_bond_capped);

        DEM_KDEM_with_damage_parallel_bond_capped() {}

        ~DEM_KDEM_with_damage_parallel_bond_capped() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        void Check(Properties::Pointer pProp) const override;

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        double GetContactSigmaMax(SphericContinuumParticle* element); // override;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseClassType)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseClassType)
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

} // namespace Kratos
#endif /* DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_CAPPED_H_INCLUDED defined */
