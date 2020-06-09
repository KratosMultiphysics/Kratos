#if !defined(DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_BILINEAR_H_INCLUDED)
#define DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_BILINEAR_H_INCLUDED

#include "DEM_KDEM_with_damage_parallel_bond_CL.h"
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_KDEM_with_damage_parallel_bond_bilinear : public DEM_KDEM_with_damage_parallel_bond {

        typedef DEM_KDEM_with_damage_parallel_bond BaseClassType;

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM_with_damage_parallel_bond_bilinear);

        DEM_KDEM_with_damage_parallel_bond_bilinear() {}

        ~DEM_KDEM_with_damage_parallel_bond_bilinear() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        void Check(Properties::Pointer pProp) const override;

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        void AdjustTauStrengthAndUpdatedMaxTauStrength(double& tau_strength, double& updated_max_tau_strength, const double internal_friction,
                                                       double contact_sigma, SphericContinuumParticle* element1, SphericContinuumParticle* element2) override;

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
#endif /* DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_BILINEAR_H_INCLUDED defined */
