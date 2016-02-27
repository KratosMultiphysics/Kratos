
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
#include "includes/define.h"
#include "spheric_particle.h"
#include "thermal_spheric_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application.h"
#include "utilities/openmp_utils.h"


//......................
//std::cout<<print...<<std::endl;

namespace Kratos
{

      ThermalSphericParticle::ThermalSphericParticle() : SphericContinuumParticle(){}

      ThermalSphericParticle::ThermalSphericParticle( IndexType NewId, GeometryType::Pointer pGeometry) : SphericContinuumParticle(NewId, pGeometry){}

      ThermalSphericParticle::ThermalSphericParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : SphericContinuumParticle(NewId, pGeometry, pProperties){}

      ThermalSphericParticle::ThermalSphericParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : SphericContinuumParticle(NewId, ThisNodes){}

      Element::Pointer ThermalSphericParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
      {
         return SphericContinuumParticle::Pointer(new ThermalSphericParticle(NewId, GetGeometry().Create( ThisNodes ), pProperties) );
      }

      /// Destructor.
      ThermalSphericParticle::~ThermalSphericParticle(){}
    
    void ThermalSphericParticle::CustomInitialize(){   
        SphericContinuumParticle::CustomInitialize();
        
        mSpecificHeat = GetProperties()[SPECIFIC_HEAT];
        mThermalConductivity = GetProperties()[THERMAL_CONDUCTIVITY]; 
        
        if (GetGeometry()[0].Coordinates()[1] > 10){ mTemperature    = 200.0; }
        else{ mTemperature    = 0.0;}
    }         
    
    double ThermalSphericParticle::GetTemperature(){return mTemperature;}

    void ThermalSphericParticle::ComputeConductiveHeatFlux(const ProcessInfo& r_process_info)
                                                       
        {                                                                                                                             
        KRATOS_TRY
                    
        /* Initializations */
        
        mConductiveHeatFlux = 0.0 ;
             
//        if (GetGeometry()[0].Coordinates()[1] < 0.02){   //0.15
//        mTemperature    = 0.0;
//            }
        
              
        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
            
            ThermalSphericParticle* neighbour_iterator = dynamic_cast<ThermalSphericParticle*>(mNeighbourElements[i]); 
            
            const double &other_radius            = neighbour_iterator->GetRadius();
            const double &other_temperature       = neighbour_iterator->GetTemperature();
            
            double rmin = GetRadius();
            if(other_radius < GetRadius()) rmin = other_radius;
            double calculation_area = KRATOS_M_PI*rmin*rmin;

            array_1d<double, 3 > other_to_me_vect = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
            
            double distance = DEM_MODULUS_3(other_to_me_vect);
            double inv_distance = 1.0 / distance;  
            mConductiveHeatFlux += - mThermalConductivity * inv_distance * calculation_area * (mTemperature - other_temperature);                                                      
        }       //for each neighbor
        if(r_process_info[TIME_STEPS]<2){KRATOS_WATCH(mConductiveHeatFlux) };
//        if (GetGeometry()[0].Coordinates()[1] > 0.29){   //0.15
//        mConductiveHeatFlux    = 1e-3;
//            }
                
        KRATOS_CATCH("")         
       
      }         //ComputeHeatFluxes  
    
    
   void ThermalSphericParticle::ComputeConvectiveHeatFlux(const ProcessInfo& r_process_info)
                                                       
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
//            ThermalSphericParticle* neighbour_iterator = dynamic_cast<ThermalSphericParticle*>(mNeighbourElements[i]); 
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
    
    void ThermalSphericParticle::CalculateRightHandSide(VectorType& r_right_hand_side_vector, 
                                                        ProcessInfo& r_current_process_info,
                                                        double dt, 
                                                        const array_1d<double,3>& gravity, int search_control)
    {
        SphericContinuumParticle::CalculateRightHandSide( r_right_hand_side_vector,  r_current_process_info, dt,  gravity, search_control);
        ThermalSphericParticle::ComputeConductiveHeatFlux(r_current_process_info);        
    }
          
    void ThermalSphericParticle::FinalizeSolutionStep(ProcessInfo& r_process_info) {
            SphericContinuumParticle::FinalizeSolutionStep(r_process_info);
            ThermalSphericParticle::UpdateTemperature(r_process_info);  
    }
    
    void ThermalSphericParticle::UpdateTemperature(const ProcessInfo& r_process_info) { 
        
            double thermal_inertia = GetMass() * mSpecificHeat;
            
            double dt = r_process_info[DELTA_TIME];
            double temperature_increment = mConductiveHeatFlux / thermal_inertia * dt;
            
            mTemperature += temperature_increment;
            GetGeometry()[0].GetSolutionStepValue(TEMPERATURE) = mTemperature; 
            GetGeometry()[0].GetSolutionStepValue(HEATFLUX) = mConductiveHeatFlux;
    }

}  // namespace Kratos.