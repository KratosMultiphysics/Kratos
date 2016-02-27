#ifndef DEM_D_DMT_COHESIVE_LAW_H
#define DEM_D_DMT_COHESIVE_LAW_H

/* Project includes */
//#include "../custom_elements/spheric_continuum_particle.h"
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    class DEM_D_DMT_Cohesive_Law : public DEMDiscontinuumConstitutiveLaw {
        
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_DMT_Cohesive_Law);

        DEM_D_DMT_Cohesive_Law();

        void Initialize(const ProcessInfo& r_process_info);

        void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

        virtual ~DEM_D_DMT_Cohesive_Law();

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const;

        double CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation = 0.0);
        double CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, DEMWall* const wall, const double indentation = 0.0);

    private:
        
        friend class Serializer;

        virtual void save(Serializer& rSerializer) const {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.save("MyMemberName", myMember);
        }

        virtual void load(Serializer& rSerializer) {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.load("MyMemberName", myMember);
        }
        
    }; //DEM_D_DMT_Cohesive_Law(const DEM_D_DMT_Cohesive_Law&);
    
} // Namespace Kratos

#endif // DEM_D_DMT_COHESIVE_LAW_H defined
