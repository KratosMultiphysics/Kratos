#ifndef DEM_D_JKR_COHESIVE_LAW_H
#define DEM_D_JKR_COHESIVE_LAW_H

/* Project includes */
//#include "../custom_elements/spheric_continuum_particle.h"
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_D_JKR_Cohesive_Law : public DEMDiscontinuumConstitutiveLaw {
        
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_JKR_Cohesive_Law);

        DEM_D_JKR_Cohesive_Law();

        void Initialize(const ProcessInfo& r_process_info) override;

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) const override;

        virtual ~DEM_D_JKR_Cohesive_Law();
        

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        double CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation) override;
        double CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, Condition* const wall, const double indentation) override;

    private:
        
        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override{
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.save("MyMemberName", myMember);
        }

        virtual void load(Serializer& rSerializer) override{
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.load("MyMemberName", myMember);
        }
    }; //DEM_D_JKR_Cohesive_Law(const DEM_D_JKR_Cohesive_Law&);
} // Namespace Kratos

#endif // DEM_D_JKR_COHESIVE_LAW_H defined
