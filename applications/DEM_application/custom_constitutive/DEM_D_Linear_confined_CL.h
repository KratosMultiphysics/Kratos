

#if !defined(DEM_D_LINEAR_CONFINED_CL_H_INCLUDED)
#define DEM_D_LINEAR_CONFINED_CL_H_INCLUDED

#include <string>
#include <iostream>
#include "DEM_D_Linear_viscous_Coulomb_CL.h"

namespace Kratos {
    
    class SphericParticle;

    class DEM_D_Linear_confined : public DEM_D_Linear_viscous_Coulomb { 
    
    public:
        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Linear_confined);

        DEM_D_Linear_confined() {}

        ~DEM_D_Linear_confined() {}

        //void Initialize(const ProcessInfo& r_process_info);

        void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;
        
        double CalculateNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation, double LocalCoordSystem[3][3]);
            
        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const;  

    }; //class DEM_D_Linear_confined

} /* namespace Kratos.*/
#endif /* DEM_D_Linear_CONFINED_CL_H_INCLUDED  defined */

