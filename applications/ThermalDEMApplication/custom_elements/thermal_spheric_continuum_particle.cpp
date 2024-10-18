//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Ferran Arrufat
//

// System includes

// External includes

// Project includes
#include "thermal_spheric_continuum_particle.h"

namespace Kratos
{
  ThermalSphericContinuumParticle::~ThermalSphericContinuumParticle() {}
  
  void ThermalSphericContinuumParticle::Initialize(const ProcessInfo& r_process_info) {
    SphericContinuumParticle::Initialize(r_process_info);

    mSpecificHeat        = GetProperties()[SPECIFIC_HEAT];
    mThermalConductivity = GetProperties()[THERMAL_CONDUCTIVITY];
  }

  const double& ThermalSphericContinuumParticle::GetTemperature() {
    return GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
  }

  void ThermalSphericContinuumParticle::SetTemperature(const double temperature) {
    GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE) = temperature;
  }

  const double& ThermalSphericContinuumParticle::GetAmbientTemperature() {
    return GetGeometry()[0].FastGetSolutionStepValue(AMBIENT_TEMPERATURE);
  }

  void ThermalSphericContinuumParticle::ComputeContactArea(const double rmin, double indentation, double& calculation_area) {
    calculation_area = Globals::Pi * rmin * rmin;
  }

  void ThermalSphericContinuumParticle::ComputeConductiveHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    mConductiveHeatFlux = 0.0;

    for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
      if (mNeighbourElements[i] == NULL) continue;
      ThermalSphericContinuumParticle* neighbour_iterator = dynamic_cast<ThermalSphericContinuumParticle*>(mNeighbourElements[i]);

      const double& other_radius      = neighbour_iterator->GetRadius();
      const double& other_temperature = neighbour_iterator->GetTemperature();

      double rmin = GetRadius();
      if (other_radius < GetRadius()) rmin = other_radius;
      
      array_1d<double, 3 > other_to_me_vect = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
      double distance     = DEM_MODULUS_3(other_to_me_vect);
      double inv_distance = 1.0 / distance;
      double radius_sum   = GetRadius() + other_radius;
      double indentation  = radius_sum - distance;

      double calculation_area = 0;
      ComputeContactArea(rmin, indentation, calculation_area);
      mConductiveHeatFlux += -mThermalConductivity * inv_distance * calculation_area * (GetTemperature() - other_temperature);
    }

    KRATOS_CATCH("")
  }

  void ThermalSphericContinuumParticle::ComputeConvectiveHeatFlux(const ProcessInfo& r_process_info) {
    //        KRATOS_TRY
    //
    //        /* Initializations */
    //        double ConvectiveHeatFlux                         = 0.0;
    //        double ambient_temperature                        = 0.0;
    //        double convective_heat_transfer_coefficient       = 5000;   // for water 5000 W/m2.K
    //
    //        if particle_is_boundary{
    //        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
    //
    //            ThermalSphericContinuumParticle<TBaseElement>* neighbour_iterator = dynamic_cast<ThermalSphericContinuumParticle<TBaseElement>*>(mNeighbourElements[i]);
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
    //            mConvectiveHeatFlux = - convective_heat_transfer_coefficient * boundary_particle_surface_area * (temperature - ambient_temperature);
    //
    //        }       //for each neighbor
    //        }
    //        KRATOS_CATCH("")
  }

  void ThermalSphericContinuumParticle::CalculateRightHandSide(const ProcessInfo& r_current_process_info, double dt, const array_1d<double, 3>& gravity) {
    SphericContinuumParticle::CalculateRightHandSide(r_current_process_info, dt, gravity);
    ComputeConductiveHeatFlux(r_current_process_info);
  }

  void ThermalSphericContinuumParticle::FinalizeSolutionStep(const ProcessInfo& r_process_info) {
    SphericContinuumParticle::FinalizeSolutionStep(r_process_info);
    UpdateTemperature(r_process_info);
    mPreviousTemperature = GetTemperature();
  }

  void ThermalSphericContinuumParticle::UpdateTemperature(const ProcessInfo& r_process_info) {
    double thermal_inertia = GetMass() * mSpecificHeat;
    double dt = r_process_info[DELTA_TIME];

    if (mSpecificHeat > 0) {// putting this condition to avoid issue, when mSpecificHeat is equal zero, which cause the wrong temperature_increment calculation
      double temperature_increment = mConductiveHeatFlux / thermal_inertia * dt;
      SetTemperature(GetTemperature() + temperature_increment);
      GetGeometry()[0].GetSolutionStepValue(HEATFLUX) = mConductiveHeatFlux;
    }
  }

  void ThermalSphericContinuumParticle::UpdateTemperatureDependentRadius(const ProcessInfo& r_process_info) {
    double thermal_alpha  = GetProperties()[THERMAL_EXPANSION_COEFFICIENT];
    double relative_temp  = GetTemperature() - GetAmbientTemperature(); // temp in Kelvin
    double updated_radius = GetRadius() * (1 + thermal_alpha * relative_temp);
    SetRadius(updated_radius);
  }

  void ThermalSphericContinuumParticle::UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion(const ProcessInfo& r_process_info, double& thermalDeltDisp, double& thermalRelVel, ThermalSphericContinuumParticle* element2) {
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

  void ThermalSphericContinuumParticle::RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(const  ProcessInfo& r_process_info,
                                                                                                       double DeltDisp[3],  //IN GLOBAL AXES
                                                                                                       double RelVel[3],    //IN GLOBAL AXES
                                                                                                       double OldLocalCoordSystem[3][3],
                                                                                                       double LocalCoordSystem[3][3],
                                                                                                       SphericParticle* neighbour_iterator) {
    //        double thermalDeltDisp = 0;
    //        double thermalRelVel   = 0;
    //        ThermalSphericContinuumParticle<TBaseElement>* thermal_neighbour_iterator = dynamic_cast<ThermalSphericContinuumParticle<TBaseElement>*>(neighbour_iterator);
    //        UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion(r_process_info, thermalDeltDisp, thermalRelVel, thermal_neighbour_iterator);
    //
    //        double LocalRelVel[3] = {0.0};
    //        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel); //TODO: can we do this in global axes directly?
    //
    //        //LocalRelVel[2] -= thermalRelVel;
    //
    //        GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalRelVel, RelVel);
  }

} // namespace Kratos.
