//
//   Project Name:                     ThermalDEM $
//   Last Modified by:    $Author: Ferran Arrufat $
//   Date:                $Date:    February 2015 $
//   Revision:            $Revision:      1.0.0.0 $
//

// System includes
#include <string>
#include <iostream>

// Project includes
#include "thermal_spheric_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{
  // Constructor/Destructor methods

  template <class TBaseElement>
  ThermalSphericParticle<TBaseElement>::~ThermalSphericParticle() {}

  // Get/Set methods

  template <class TBaseElement>
  const double& ThermalSphericParticle<TBaseElement>::GetParticleTemperature() {
    return GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
  }

  template <class TBaseElement>
  const double& ThermalSphericParticle<TBaseElement>::GetInterstitialFluidTemperature() {
    return GetGeometry()[0].FastGetSolutionStepValue(AMBIENT_TEMPERATURE);
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::SetParticleTemperature(const double temperature) {
    GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE) = temperature;
  }

  // Initialization methods

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::Initialize(const ProcessInfo& r_process_info) {
    TBaseElement::Initialize(r_process_info);
    mSpecificHeat        = GetProperties()[SPECIFIC_HEAT];
    mThermalConductivity = GetProperties()[THERMAL_CONDUCTIVITY];
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::InitializeSolutionStep(const ProcessInfo& r_process_info) {
    TBaseElement::InitializeSolutionStep(r_process_info);
  }

  // Calculate right hand side (forces and heat fluxes)

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::CalculateRightHandSide(const ProcessInfo& r_current_process_info, double dt, const array_1d<double, 3>& gravity) {
    // Force components
    TBaseElement::CalculateRightHandSide(r_current_process_info, dt, gravity);
    
    // Heat flux components
    ComputeHeatFluxes(r_current_process_info);
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeHeatFluxes(const ProcessInfo& r_process_info) {
    // Initialize heat flux contributions
    mConductiveHeatFlux = 0.0;
    mConvectiveHeatFlux = 0.0;
    mTotalHeatFlux      = 0.0;

    // Direct conduction
    ComputeBallToBallDirectConductionHeatFlux(r_process_info);
    ComputeBallToRigidFaceDirectConductionHeatFlux(r_process_info);

    // Convection
    ComputeConvectiveHeatFlux(r_process_info);

    // Other sources

    // Sum up contributions
    mTotalHeatFlux = mConductiveHeatFlux + mConvectiveHeatFlux;
  }

  // Compute heat fluxes components

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeBallToBallDirectConductionHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Loop over neighbour elements
    for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
      if (mNeighbourElements[i] == NULL) continue;
      ThermalSphericParticle<TBaseElement>* neighbour_iterator = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(mNeighbourElements[i]);

      // Get particle properties
      double this_temp   = GetParticleTemperature();
      double this_radius = GetRadius();

      // Get neighbour properties
      double other_temp   = neighbour_iterator->GetParticleTemperature();
      double other_radius = neighbour_iterator->GetRadius();

      // Compute minimum radius
      double rmin = std::min(this_radius,other_radius);

      // Compute interaction parameters
      array_1d<double,3> other_to_this_vect = GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
      double distance     = DEM_MODULUS_3(other_to_this_vect);
      double inv_distance = 1.0 / distance;
      double indentation  = this_radius + other_radius - distance;

      // Compute contact area
      double contact_area = 0.0;
      ComputeContactArea(rmin, indentation, contact_area);

      // Compute heat flux
      mConductiveHeatFlux += mThermalConductivity * inv_distance * contact_area * (other_temp - this_temp);
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeBallToRigidFaceDirectConductionHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    std::vector<DEMWall*>& rNeighbours = mNeighbourRigidFaces;

    // Loop over neighbour walls
    for (unsigned int i = 0; i < rNeighbours.size(); i++) {
      DEMWall* wall = rNeighbours[i];
      if (wall == NULL) continue;
      
      // Get particle properties
      double particle_temp = GetParticleTemperature();

      // Get wall temperature as the average of its nodes
      double wall_temp = 0.0;
      double n_nodes = wall->GetGeometry().size();
      for (unsigned int i = 0; i < n_nodes; i++) {
        double node_temp = wall->GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE);
        wall_temp += node_temp;
      }
      wall_temp /= n_nodes;

      // Compute contact area
      double contact_radius = 0.1;

      // Compute heat flux (Batchelor & OBrien model)
      mConductiveHeatFlux += 2 * mThermalConductivity * contact_radius * (wall_temp - particle_temp);
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeConvectiveHeatFlux(const ProcessInfo& r_process_info) {
    //KRATOS_TRY

    ///* Initializations */

    //double ConvectiveHeatFlux                   = 0.0;
    //double ambient_temperature                  = 0.0;
    //double convective_heat_transfer_coefficient = 5000;  // for water 5000 W/m2.K

    //if particle_is_boundary {
    //  for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {

    //    ThermalSphericParticle<TBaseElement>* neighbour_iterator = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(mNeighbourElements[i]);

    //    array_1d<double, 3 > other_to_me_vect = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();

    //    double distance = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
    //                      other_to_me_vect[1] * other_to_me_vect[1] +
    //                      other_to_me_vect[2] * other_to_me_vect[2]);

    //    double ThermalConductivity = 50;
    //    double inv_distance = 1/distance;
    //    //double inv_distance = 1/(GetRadius()+other_radius);

    //    mConvectiveHeatFlux = - convective_heat_transfer_coefficient * boundary_particle_surface_area * (temperature - ambient_temperature);

    //  }       //for each neighbor
    //}
    //KRATOS_CATCH("")
  }

  // Auxiliary computation methods

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeContactArea(const double rmin, double indentation, double& calculation_area) {
    calculation_area = Globals::Pi*rmin*rmin;
  }

  // Update methods

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::UpdateTemperature(const ProcessInfo& r_process_info) {
    // Condition to avoid issue when mSpecificHeat is equal zero, which cause the wrong temperature_increment calculation
    if (mSpecificHeat > 0) {
      // Compute new temperature
      double dt              = r_process_info[DELTA_TIME];
      double thermal_inertia = GetMass() * mSpecificHeat;
      double temp_increment  = mTotalHeatFlux / thermal_inertia * dt;
      double temp_new        = GetParticleTemperature() + temp_increment;
      
      // Set new temperature
      SetParticleTemperature(temp_new);
    }
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::UpdateTemperatureDependentRadius(const ProcessInfo& r_process_info) {
    double this_temp      = GetParticleTemperature();
    double fluid_temp     = GetInterstitialFluidTemperature();
    double relative_temp  = this_temp - fluid_temp; // temp in Kelvin
    double thermal_alpha  = GetProperties()[THERMAL_EXPANSION_COEFFICIENT];
    double updated_radius = GetRadius() * (1 + thermal_alpha * relative_temp);
    SetRadius(updated_radius);
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion(const ProcessInfo& r_process_info, double& thermalDeltDisp, double& thermalRelVel, ThermalSphericParticle<TBaseElement>* element2) {
    //thermalRelVel = 0;
    //double temperature = GetParticleTemperature();
    //double other_temperature = element2->GetParticleTemperature();
    //double previous_temperature = mPreviousTemperature;
    //double previous_other_temperature = element2->mPreviousTemperature;
    //double thermal_alpha = GetProperties()[THERMAL_EXPANSION_COEFFICIENT];
    //double updated_radius = GetRadius();
    //double updated_other_radius = element2->GetRadius();
    //double dt = r_process_info[DELTA_TIME];
    //double temperature_increment_elem1 = temperature - previous_temperature;
    //double temperature_increment_elem2 = other_temperature - previous_other_temperature;
    //thermalDeltDisp = updated_radius * thermal_alpha * temperature_increment_elem1;
    //thermalDeltDisp = thermalDeltDisp + updated_other_radius * thermal_alpha * temperature_increment_elem2;
    //thermalRelVel = updated_radius * thermal_alpha * temperature_increment_elem1 / dt;
    //thermalRelVel = thermalRelVel + updated_other_radius * thermal_alpha * temperature_increment_elem2 / dt;
  }

  // Finalization methods

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::FinalizeSolutionStep(const ProcessInfo& r_process_info) {
    TBaseElement::FinalizeSolutionStep(r_process_info);
    UpdateTemperature(r_process_info);
    mPreviousTemperature = GetParticleTemperature();
    GetGeometry()[0].GetSolutionStepValue(HEATFLUX) = mTotalHeatFlux;
  }

  // Others

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(const ProcessInfo& r_process_info,
                                                                                                            double DeltDisp[3], //IN GLOBAL AXES
                                                                                                            double RelVel[3], //IN GLOBAL AXES
                                                                                                            double OldLocalCoordSystem[3][3],
                                                                                                            double LocalCoordSystem[3][3],
                                                                                                            SphericParticle* neighbour_iterator) {
    //double thermalDeltDisp = 0;
    //double thermalRelVel = 0;
    //ThermalSphericParticle<TBaseElement>* thermal_neighbour_iterator = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(neighbour_iterator);
    //UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion(r_process_info, thermalDeltDisp, thermalRelVel, thermal_neighbour_iterator);
    //double LocalRelVel[3] = {0.0};
    //GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel); //TODO: can we do this in global axes directly?
    ////LocalRelVel[2] -= thermalRelVel;
    //GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalRelVel, RelVel);
  }

  //Explicit Instantiation
  template class ThermalSphericParticle<SphericParticle>;
  template class ThermalSphericParticle<SphericContinuumParticle>;

} // namespace Kratos
