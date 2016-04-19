
#include <iostream>

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
                                                          double ElasticLocalRotationalMoment[3],
                                                          double ViscoLocalRotationalMoment[3]) {
        KRATOS_TRY
        
        double fabric_coefficient = element->GetProperties()[FABRIC_COEFFICIENT];
        
        DEM_KDEM::ComputeParticleRotationalMoments(element, neighbor, equiv_young, distance, calculation_area, LocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment);
        
        DEM_MULTIPLY_BY_SCALAR_3(ElasticLocalRotationalMoment, fabric_coefficient);
        DEM_MULTIPLY_BY_SCALAR_3(ViscoLocalRotationalMoment, fabric_coefficient);
        
        KRATOS_CATCH("")
    }
    
    
    void DEM_KDEMFabric::AddPoissonContribution(const double equiv_poisson, double LocalCoordSystem[3][3], double& normal_force, 
                                          double calculation_area, Matrix* mSymmStressTensor, SphericParticle* element1, SphericParticle* element2) {
        KRATOS_TRY         
        KRATOS_CATCH("")
    }

} // namespace Kratos
