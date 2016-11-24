// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_Mohr_Coulomb_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_Mohr_Coulomb::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_Mohr_Coulomb(*this));
        return p_clone;
    }

    void DEM_KDEM_Mohr_Coulomb::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        std::cout << "\nAssigning DEM_KDEM_Mohr_Coulomb to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }    
    
    double DEM_KDEM_Mohr_Coulomb::LocalMaxSearchDistance(const int i, SphericContinuumParticle* element1, SphericContinuumParticle* element2) {

        Properties& element1_props = element1->GetProperties();
        Properties& element2_props = element2->GetProperties();
        const double mohr_coulomb_c = 1e6 * 0.5*(element1_props[INTERNAL_COHESION] + element2_props[INTERNAL_COHESION]);

        // calculation of equivalent young modulus
        double myYoung = element1->GetYoung();
        double other_young = element2->GetYoung();
        double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);

        const double my_radius = element1->GetRadius();
        const double other_radius = element2->GetRadius();
        double calculation_area = 0;
        Vector& vector_of_contact_areas = element1->GetValue(NEIGHBOURS_CONTACT_AREAS);

        GetContactArea(my_radius, other_radius, vector_of_contact_areas, i, calculation_area); 

        const double radius_sum = my_radius + other_radius;
        const double initial_delta = element1->GetInitialDelta(i);
        const double initial_dist = radius_sum - initial_delta;
        const double kn_el = equiv_young * calculation_area / initial_dist;        

        const double max_normal_force = mohr_coulomb_c * calculation_area;
        const double u1 = max_normal_force / kn_el;
        return u1;
    }
    
    void DEM_KDEM_Mohr_Coulomb::CheckFailure(const int i_neighbour_count, SphericContinuumParticle* element1, SphericContinuumParticle* element2){
        
        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
        
        if (failure_type == 0) {  
            Matrix average_stress_tensor = ZeroMatrix(3,3);       
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    average_stress_tensor(i,j) = 0.5 * ((*(element1->mSymmStressTensor))(i,j) + (*(element2->mSymmStressTensor))(i,j));
                }
            }

            Vector principal_stresses(3);
            noalias(principal_stresses) = AuxiliaryFunctions::EigenValuesDirectMethod(average_stress_tensor);
            
            Properties& element1_props = element1->GetProperties();
            Properties& element2_props = element2->GetProperties();
            
            const double mohr_coulomb_c = 1e6 * 0.5*(element1_props[INTERNAL_COHESION] + element2_props[INTERNAL_COHESION]);
            const double mohr_coulomb_phi = 0.5 * (element1_props[INTERNAL_FRICTION_ANGLE] + element2_props[INTERNAL_FRICTION_ANGLE]);
            const double mohr_coulomb_phi_in_radians = mohr_coulomb_phi * KRATOS_M_PI / 180.0; 
            const double sinphi = sin(mohr_coulomb_phi_in_radians);
            const double cosphi = cos(mohr_coulomb_phi_in_radians);
            
            const double function_value = (principal_stresses[0] - principal_stresses[2]) + (principal_stresses[0] + principal_stresses[2]) * sinphi - 2.0 * mohr_coulomb_c * cosphi;
            
            if(function_value > 0) {                
                failure_type = 4;
            }
                         
        }    
        
    }
     

} // namespace Kratos
