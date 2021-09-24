#if !defined(DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_HERTZ_H_INCLUDED)
#define DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_HERTZ_H_INCLUDED

/* Project includes */
#include "DEM_KDEM_with_damage_parallel_bond_CL.h"
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_KDEM_with_damage_parallel_bond_Hertz : public DEM_KDEM_with_damage_parallel_bond {

        typedef DEM_KDEM_with_damage_parallel_bond BaseClassType;

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM_with_damage_parallel_bond_Hertz);

        DEM_KDEM_with_damage_parallel_bond_Hertz() {}

        ~DEM_KDEM_with_damage_parallel_bond_Hertz() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        void SetConstitutiveLawInPropertiesWithParameters(Properties::Pointer pProp, const Parameters& parameters, bool verbose) override;

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        void CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double bonded_equiv_young,
                                             double equiv_poisson, double calculation_area,
                                             SphericContinuumParticle* element1, SphericContinuumParticle* element2, double indentation) override;

        void ComputeNormalUnbondedForce(double indentation) override;

        double LocalMaxSearchDistance(const int i, SphericContinuumParticle* element1, SphericContinuumParticle* element2) override;

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
#endif /* DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_HERTZ_H_INCLUDED defined */
