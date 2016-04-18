// System includes
#include <iostream>

// Project includes
#include "DEM_application.h"
#include "DEM_KDEM_fabric_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEMFabric::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEMFabric(*this));
        return p_clone;
    }


    void DEM_KDEMFabric::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        std::cout << "\nAssigning DEM_KDEMFabric to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }    
                                  
    
    void DEM_KDEMFabric::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                    SphericContinuumParticle* neighbor,
                                                    double equiv_young,
                                                    double distance,
                                                    double calculation_area,
                                                    double LocalCoordSystem[3][3],
                                                    array_1d<double, 3>& mContactMoment) {
        KRATOS_TRY
        
        DEM_KDEM::ComputeParticleRotationalMoments(element, neighbor, equiv_young, distance, calculation_area, LocalCoordSystem, mContactMoment);
        mContactMoment *= 0.1; // TODO: Hardcoded the reduction of a 90% of the flection in relation to KDEM
        
        KRATOS_CATCH("")
    }
    
    
    void DEM_KDEMFabric::AddPoissonContribution(const double equiv_poisson, double LocalCoordSystem[3][3], double& normal_force, 
                                          double calculation_area, Matrix* mSymmStressTensor, SphericParticle* element1, SphericParticle* element2) {
        KRATOS_TRY         
        KRATOS_CATCH("")
    }

} // namespace Kratos
