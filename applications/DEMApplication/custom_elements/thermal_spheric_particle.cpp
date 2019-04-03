
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
#include "utilities/openmp_utils.h"


//......................
//KRATOS_INFO("DEM") <<print...<<std::endl;

namespace Kratos
{
    template< class TBaseElement >
    ThermalSphericParticle<TBaseElement>::~ThermalSphericParticle() {}


    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::Initialize(const ProcessInfo& r_process_info){
        TBaseElement::Initialize(r_process_info);
         // We can put here UpdateTemperatureDependentRadius()
        mSpecificHeat = GetProperties()[SPECIFIC_HEAT];
        mThermalConductivity = GetProperties()[THERMAL_CONDUCTIVITY];
        //mPreviousTemperature = 250;	// initial temperature in the first time step
        //SetTemperature(mPreviousTemperature);// initial temperature in the first time step
        if (GetGeometry()[0].Coordinates()[1] > 2){ GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE)    = 200.0; }
        else{ GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE)    = 0.0;}
    }

    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::InitializeSolutionStep(ProcessInfo& r_process_info)
    {
            TBaseElement::InitializeSolutionStep(r_process_info);
            //UpdateTemperatureDependentRadius(r_process_info);
    }

    template< class TBaseElement >
    const double& ThermalSphericParticle<TBaseElement>::GetTemperature(){return GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);}

    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::SetTemperature(const double temperature){GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE)=temperature;}


    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::ComputeContactArea(const double rmin, double indentation, double& calculation_area) {
            calculation_area = Globals::Pi*rmin*rmin;
    }

    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::ComputeConductiveHeatFlux(const ProcessInfo& r_process_info) {
        KRATOS_TRY
        mConductiveHeatFlux = 0.0 ;

//        if (GetGeometry()[0].Coordinates()[1] < 0.02){   //0.15
//        mTemperature    = 0.0;
//            }

        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
            if(mNeighbourElements[i] == NULL) continue;

            ThermalSphericParticle<TBaseElement>* neighbour_iterator = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(mNeighbourElements[i]);
            const double &other_radius = neighbour_iterator->GetRadius();
            const double &other_temperature = neighbour_iterator->GetTemperature();

            double rmin = GetRadius();
            if (other_radius < GetRadius()) rmin = other_radius;
            double calculation_area = 0;
            array_1d<double, 3 > other_to_me_vect = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();

            double distance = DEM_MODULUS_3(other_to_me_vect);
            double inv_distance = 1.0 / distance;
            double radius_sum = GetRadius() + other_radius;
            double indentation = radius_sum - distance;

            ComputeContactArea(rmin, indentation, calculation_area);
            mConductiveHeatFlux += -mThermalConductivity * inv_distance * calculation_area * (GetTemperature() - other_temperature);
        }       //for each neighbor

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
                                                        const array_1d<double,3>& gravity, int search_control) {
        TBaseElement::CalculateRightHandSide(r_current_process_info, dt,  gravity, search_control);
        ComputeConductiveHeatFlux(r_current_process_info);
    }


    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::FinalizeSolutionStep(ProcessInfo& r_process_info) {
            TBaseElement::FinalizeSolutionStep(r_process_info);
            UpdateTemperature(r_process_info);
            mPreviousTemperature = GetTemperature();
    }

    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::UpdateTemperature(const ProcessInfo& r_process_info) {
        double thermal_inertia = GetMass() * mSpecificHeat;
        double dt = r_process_info[DELTA_TIME];
        if (mSpecificHeat>0) {// putting this condition to avoid issue, when mSpecificHeat is equal zero, which cause the wrong temperature_increment calculation
            double temperature_increment = mConductiveHeatFlux / thermal_inertia * dt;
            SetTemperature(GetTemperature() + temperature_increment);
            GetGeometry()[0].GetSolutionStepValue(HEATFLUX) = mConductiveHeatFlux;
        }
    }

    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::UpdateTemperatureDependentRadius(const ProcessInfo& r_process_info)
    {
        double temperature = GetTemperature();
        double radius = GetRadius();
        double thermal_alpha = GetProperties()[THERMAL_EXPANSION_COEFFICIENT];
        const int RT = 293; // room temperature // is there already a variable for RT?
        double relative_temp = temperature - RT; // temp in Kelvin
        double updated_radius = radius * (1 + thermal_alpha*relative_temp);
        SetRadius(updated_radius);
    }

    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion(const ProcessInfo& r_process_info, double& thermalDeltDisp, double& thermalRelVel, ThermalSphericParticle<TBaseElement>* element2)
    {
//        thermalRelVel = 0;
//        double temperature = GetTemperature();
//        double other_temperature = element2->GetTemperature();
//        double previous_temperature = mPreviousTemperature;
//        double previous_other_temperature = element2->mPreviousTemperature;
//        double thermal_alpha = GetProperties()[THERMAL_EXPANSION_COEFFICIENT];
//        double updated_radius = GetRadius();
//        double updated_other_radius = element2->GetRadius();

//        double dt = r_process_info[DELTA_TIME];
//        double temperature_increment_elem1 = temperature - previous_temperature;
//        double temperature_increment_elem2 = other_temperature - previous_other_temperature;
//        thermalDeltDisp = updated_radius * thermal_alpha * temperature_increment_elem1;
//        thermalDeltDisp = thermalDeltDisp + updated_other_radius * thermal_alpha * temperature_increment_elem2;
//        thermalRelVel = updated_radius * thermal_alpha * temperature_increment_elem1 / dt;
//        thermalRelVel = thermalRelVel + updated_other_radius * thermal_alpha * temperature_increment_elem2 / dt;
    }

    template< class TBaseElement >
    void ThermalSphericParticle<TBaseElement>::RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(const ProcessInfo& r_process_info,
                                                                                                            double DeltDisp[3], //IN GLOBAL AXES
                                                                                                            double RelVel[3], //IN GLOBAL AXES
                                                                                                            double OldLocalCoordSystem[3][3],
                                                                                                            double LocalCoordSystem[3][3],
                                                                                                            SphericParticle* neighbour_iterator) {
//        double thermalDeltDisp = 0;
//        double thermalRelVel = 0;
//        ThermalSphericParticle<TBaseElement>* thermal_neighbour_iterator = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(neighbour_iterator);
//        UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion(r_process_info, thermalDeltDisp, thermalRelVel, thermal_neighbour_iterator);

//        double LocalRelVel[3]          = {0.0};
//        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel); //TODO: can we do this in global axes directly?

//        //LocalRelVel[2] -= thermalRelVel;

//        GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalRelVel, RelVel);
    }

    template class ThermalSphericParticle<SphericParticle>; //Explicit Instantiation
    template class ThermalSphericParticle<SphericContinuumParticle>; //Explicit Instantiation
}  // namespace Kratos.
