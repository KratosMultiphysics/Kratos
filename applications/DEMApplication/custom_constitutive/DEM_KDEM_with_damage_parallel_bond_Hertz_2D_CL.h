#if !defined(DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_HERTZ_2D_H_INCLUDED)
#define DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_HERTZ_2D_H_INCLUDED

/* Project includes */
#include "DEM_KDEM_with_damage_parallel_bond_Hertz_CL.h"
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_KDEM_with_damage_parallel_bond_Hertz_2D : public DEM_KDEM_with_damage_parallel_bond_Hertz {

        typedef DEM_KDEM_with_damage_parallel_bond_Hertz BaseClassType;

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM_with_damage_parallel_bond_Hertz_2D);

        DEM_KDEM_with_damage_parallel_bond_Hertz_2D() {}

        ~DEM_KDEM_with_damage_parallel_bond_Hertz_2D() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        void SetConstitutiveLawInPropertiesWithParameters(Properties::Pointer pProp, const Parameters& parameters, bool verbose) override;

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        void CalculateContactArea(double radius, double other_radius, double& calculation_area) override;

        void CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young,
                                       double equiv_poisson, double calculation_area, SphericContinuumParticle* element1,
                                       SphericContinuumParticle* element2, double indentation) override;
        
        void ComputeNormalUnbondedForce(double indentation) override;

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
#endif /* DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_HERTZ_2D_H_INCLUDED defined */
