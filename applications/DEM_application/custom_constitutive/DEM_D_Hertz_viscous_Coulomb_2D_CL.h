//Authors: M.A. Celigueta and S. Latorre (CIMNE)
//   Date: July 2015

#if !defined(DEM_D_HERTZ_VISCOUS_COULOMB_2D_CL_H_INCLUDED)
#define  DEM_D_HERTZ_VISCOUS_COULOMB_2D_CL_H_INCLUDED

#include <string>
#include <iostream>

#include "DEM_application.h"
#include "includes/define.h"
#include "DEM_D_Hertz_viscous_Coulomb_CL.h"

namespace Kratos {
    
    class SphericParticle;

    class DEM_D_Hertz_viscous_Coulomb2D : public DEM_D_Hertz_viscous_Coulomb {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Hertz_viscous_Coulomb2D);

        DEM_D_Hertz_viscous_Coulomb2D() {
        }
      
        void Initialize(const ProcessInfo& r_process_info);         

        void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

        ~DEM_D_Hertz_viscous_Coulomb2D() {
        }

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const;      
        
        void InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation);  
        
        void InitializeContactWithFEM(SphericParticle* const element, DEMWall* const wall, const double indentation, const double ini_delta = 0.0);
                
                            
    private:

        friend class Serializer;
        
        virtual void save(Serializer& rSerializer) const {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

} /* namespace Kratos.*/
#endif /* DEM_D_HERTZ_VISCOUS_COULOMB_2D_CL_H_INCLUDED  defined */
