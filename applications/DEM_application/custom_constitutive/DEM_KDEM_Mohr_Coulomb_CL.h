
#if !defined(DEM_KDEM_MOHR_COULOMB_CL_H_INCLUDED)
#define  DEM_KDEM_MOHR_COULOMB_CL_H_INCLUDED

/* Project includes */
#include "DEM_continuum_constitutive_law.h"
#include "DEM_KDEM_Rankine_CL.h"


namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_KDEM_Mohr_Coulomb : public DEM_KDEM_Rankine {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM_Mohr_Coulomb);

        DEM_KDEM_Mohr_Coulomb() {
        }    

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) const override;

        ~DEM_KDEM_Mohr_Coulomb() {
        }

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;
        
        double LocalMaxSearchDistance(const int i, SphericContinuumParticle* element1, SphericContinuumParticle* element2) override;
        void CheckFailure(const int i_neighbour_count, SphericContinuumParticle* element1, SphericContinuumParticle* element2) override;
        bool CheckRequirementsOfStressTensor() override;            

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override{
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEM_KDEM)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override{
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEM_KDEM)
                    //rSerializer.load("MyMemberName",myMember);
        }

    };

} /* namespace Kratos.*/
#endif /* DEM_KDEM_MOHR_COULOMB_H_INCLUDED  defined */
