
//   Project Name:        Kratos
//   Last Modified by:    $Author: Ferran Arrufat $
//   Date:                $Date: 2016-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes
#include <string>
#include <iostream> 

// External includes


// Project includes
#include "thermal_spheric_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application.h"
#include "utilities/openmp_utils.h"


//......................
//std::cout<<print...<<std::endl;

namespace Kratos
{      
    
    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::CustomInitialize(){   
        TBaseElement::CustomInitialize();
         // We can put here UpdateTemperatureDependentRadius()
        mSpecificHeat = GetProperties()[SPECIFIC_HEAT];
        mThermalConductivity = GetProperties()[THERMAL_CONDUCTIVITY]; 
        if (GetGeometry()[0].Coordinates()[1] > 10){ GetTemperature() = 200.0; }
        else{ GetTemperature() = 0.0;}
    }         
    
    template< class TBaseElement >
    double& ThermalSphericParticle<TBaseElement>::GetTemperature(){return GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);}
	template< class TBaseElement >
	double& ThermalSphericParticle<TBaseElement>::GetPreviousTemperature(){ return mPrevTemperature; }
    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::ComputeConductiveHeatFlux(const ProcessInfo& r_process_info)
                                                       
        {                                                                                                                             
        KRATOS_TRY
                    
        /* Initializations */
        
        mConductiveHeatFlux = 0.0 ;
             
//        if (GetGeometry()[0].Coordinates()[1] < 0.02){   //0.15
//        mTemperature    = 0.0;
//            }
        
              
        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
            
            ThermalSphericParticle<TBaseElement>* neighbour_iterator = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(mNeighbourElements[i]); 
            
            const double &other_radius            = neighbour_iterator->GetRadius();
            const double &other_temperature       = neighbour_iterator->GetTemperature();
            
            double rmin = GetRadius();
            if(other_radius < GetRadius()) rmin = other_radius;
            double calculation_area = KRATOS_M_PI*rmin*rmin;

            array_1d<double, 3 > other_to_me_vect = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
            
            double distance = DEM_MODULUS_3(other_to_me_vect);
            double inv_distance = 1.0 / distance;  
            mConductiveHeatFlux += - mThermalConductivity * inv_distance * calculation_area * (GetTemperature() - other_temperature);
        }       //for each neighbor
//        if (GetGeometry()[0].Coordinates()[1] > 0.29){   //0.15
//        mConductiveHeatFlux    = 1e-3;
//            }
                
        KRATOS_CATCH("")         
       
      }         //ComputeHeatFluxes  
    
    
   template< class TBaseElement > 
   void ThermalSphericParticle<TBaseElement>::ComputeConvectiveHeatFlux(const ProcessInfo& r_process_info)
                                                       
        {                                                                                                                             
//        KRATOS_TRY
//                    
//        /* Initializations */
//        
//        double ConvectiveHeatFlux                         = 0.0;
//        double ambient_temperature                        = 0.0;
//        double convective_heat_transfer_coefficient       = 5000;   // for water 5000 W/m2.K
//       
//        if particle_is_boundary{      
//        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
//            
//            ThermalSphericParticle<TBaseElement>* neighbour_iterator = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(mNeighbourElements[i]); 
//            
//            array_1d<double, 3 > other_to_me_vect = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
//            
//            double distance = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
//                    other_to_me_vect[1] * other_to_me_vect[1] +
//                    other_to_me_vect[2] * other_to_me_vect[2]);
//            
//            double ThermalConductivity = 50;
//            double inv_distance = 1/distance;
//            //double inv_distance = 1/(GetRadius()+other_radius);
//            
//            mConvectiveHeatFlux = - convective_heat_transfer_coefficient * boundary_particle_surface_area * (mTemperature - ambient_temperature);
//                                  
//        }       //for each neighbor
//        }
//        KRATOS_CATCH("")         
       
      }         //ComputeHeatFluxes      
    
    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::CalculateRightHandSide(ProcessInfo& r_current_process_info,
                                                        double dt, 
                                                        const array_1d<double,3>& gravity, int search_control)
    {
        TBaseElement::CalculateRightHandSide(r_current_process_info, dt,  gravity, search_control);
        ComputeConductiveHeatFlux(r_current_process_info);        
    }
        
    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::FinalizeSolutionStep(ProcessInfo& r_process_info) {
            TBaseElement::FinalizeSolutionStep(r_process_info);
            UpdateTemperature(r_process_info);  
			//SetValue(TEMPERATURE_OLD_IT, GetTemperature());
    }
    
    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::UpdateTemperature(const ProcessInfo& r_process_info) { 
        
            double thermal_inertia = GetMass() * mSpecificHeat;
            
            double dt = r_process_info[DELTA_TIME];
            double temperature_increment = mConductiveHeatFlux / thermal_inertia * dt;
            
			GetTemperature() += temperature_increment;
            //GetGeometry()[0].GetSolutionStepValue(TEMPERATURE) = GetTemperature();
            GetGeometry()[0].GetSolutionStepValue(HEATFLUX) = mConductiveHeatFlux;
    }

	template< class TBaseElement >
	void ThermalSphericParticle<TBaseElement>::UpdateTemperatureDependentRadius(const ProcessInfo& r_process_info)
	{
		double temperature = GetTemperature();
		//KRATOS_WATCH(temperature);
		double radius = GetRadius();
		//KRATOS_WATCH(GetRadius());
		double thermal_alpha = GetProperties()[THERMAL_EXPANSION_COEFFICIENT];
		//KRATOS_WATCH(thermal_alpha);
		const int RT = 293; // initial temperature // is there already a variable for RT?
		double relative_temp = temperature - RT; // temp in Kelvin
		//KRATOS_WATCH(relative_temp);
		double updated_radius = radius * (1 + thermal_alpha*relative_temp);
		SetRadius(updated_radius);
		//KRATOS_WATCH(GetRadius());
	}

	template< class TBaseElement >
	void ThermalSphericParticle<TBaseElement>::UpdateTemperatureDependentRadii(const ProcessInfo& r_process_info, ThermalSphericParticle<TBaseElement>* element2)
	{
		double temperature = GetTemperature();
		double other_temperature = element2->GetTemperature();
		double my_radius = GetRadius();
	    double other_radius = element2->GetRadius();
		double thermal_alpha = GetProperties()[THERMAL_EXPANSION_COEFFICIENT];
		const int RT = 293; // initial temperature // is there already a variable for RT?
		double relative_temp_1 = temperature - RT; // temp in Kelvin
		double relative_temp_2 = other_temperature - RT; // temp in Kelvin
		double updated_radius = my_radius * (1 + thermal_alpha*relative_temp_1);
		double updated_other_radius = other_radius * (1 + thermal_alpha*relative_temp_2);
		SetRadius(updated_radius);
		element2->SetRadius(updated_other_radius);
		//KRATOS_WATCH(GetRadius());
	}

	template< class TBaseElement >
	void ThermalSphericParticle<TBaseElement>::UpdateNormalRelativeVelocityDueToThermalExpansion(const ProcessInfo& r_process_info, double& thermalRelVel, ThermalSphericParticle<TBaseElement>* element2)
	{
		double temperature = GetTemperature();
		double other_temperature = element2->GetTemperature();
		double previous_temperature = GetValue(TEMPERATURE_OLD_IT);
		double previous_other_temperature = element2->GetValue(TEMPERATURE_OLD_IT);
		const int RT = 293; // initial temperature // is there already a variable for RT?
		double relative_temp_1 = temperature - RT; // temp in Kelvin
		double relative_temp_2 = other_temperature - RT; // temp in Kelvin

		double my_radius = GetRadius();
		double other_radius = element2->GetRadius();
		double thermal_alpha = GetProperties()[THERMAL_EXPANSION_COEFFICIENT];

		double updated_radius = my_radius * (1 + thermal_alpha*relative_temp_1);
		double updated_other_radius = other_radius * (1 + thermal_alpha*relative_temp_2);

		double dt = r_process_info[DELTA_TIME];
		double temperature_increment_elem1 = temperature - previous_temperature;
		double temperature_increment_elem2 = other_temperature - previous_other_temperature;
		thermalRelVel = updated_radius * thermal_alpha * temperature_increment_elem1 / dt;
		thermalRelVel = thermalRelVel + updated_other_radius * thermal_alpha * temperature_increment_elem2 / dt;
		SetValue(TEMPERATURE_OLD_IT, temperature);
	}
        
        template< class TBaseElement >
        void ThermalSphericParticle<TBaseElement>::RelativeDisplacementAndRotationOfContactPointDueToOtherReasons(const ProcessInfo& r_process_info,
                                                                                                                double DeltDisp[3], //IN GLOBAL AXES
                                                                                                                double RelVel[3], //IN GLOBAL AXES
                                                                                                                double OldLocalCoordSystem[3][3],
                                                                                                                double LocalCoordSystem[3][3], 
                                                                                                                SphericParticle* neighbour_iterator) {
            
            double thermalRelVel = 0;
            //UpdateTemperatureDependentRadii(r_process_info, neighbour_iterator);
            ThermalSphericParticle<TBaseElement>* thermal_neighbour_iterator = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(neighbour_iterator);
            UpdateNormalRelativeVelocityDueToThermalExpansion(r_process_info, thermalRelVel, thermal_neighbour_iterator);
            
            double LocalRelVel[3]          = {0.0};
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel); //TODO: can we do this in global axes directly?
            LocalRelVel[2] -= thermalRelVel;
            GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalRelVel, RelVel);                
        }

    template class ThermalSphericParticle<SphericParticle>; //Explicit Instantiation
    template class ThermalSphericParticle<SphericContinuumParticle>; //Explicit Instantiation
}  // namespace Kratos.
