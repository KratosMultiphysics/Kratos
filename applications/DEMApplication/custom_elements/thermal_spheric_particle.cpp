//
//   Project Name:                     ThermalDEM $
//   Last Modified by:    $Author: Ferran Arrufat $
//   Date:                $Date:    February 2015 $
//   Revision:            $Revision:      1.0.0.0 $
//

// System includes
#include <string>
#include <iostream>
#include <limits>

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

  //=====================================================================================================================================================================================
  // Initialization methods

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::Initialize(const ProcessInfo& r_process_info) {
    // Initialize base class
    TBaseElement::Initialize(r_process_info);

    // Set thermal flags
    this->Set(DEMFlags::HAS_DIRECT_CONDUCTION,   r_process_info[DIRECT_CONDUCTION_OPTION]);
    this->Set(DEMFlags::HAS_INDIRECT_CONDUCTION, r_process_info[INDIRECT_CONDUCTION_OPTION]);
    this->Set(DEMFlags::HAS_CONVECTION,          r_process_info[CONVECTION_OPTION]);
    this->Set(DEMFlags::HAS_RADIATION,           r_process_info[CONVECTION_OPTION]);

    // Initialize prescribed heat flux
    SetParticlePrescribedHeatFlux(0.0);
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::InitializeSolutionStep(const ProcessInfo& r_process_info) {
    // Initialize base class
    TBaseElement::InitializeSolutionStep(r_process_info);
  }

  //=====================================================================================================================================================================================
  // Calculate right hand side (forces and heat fluxes)

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::CalculateRightHandSide(const ProcessInfo& r_process_info, double dt, const array_1d<double, 3>& gravity) {
    // Force components
    TBaseElement::CalculateRightHandSide(r_process_info, dt, gravity);
    
    // Heat flux components
    ComputeHeatFluxes(r_process_info);
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeHeatFluxes(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Initialize heat fluxes contributions
    mConductiveHeatFlux = 0.0;
    mConvectiveHeatFlux = 0.0;
    mRadiativeHeatFlux  = 0.0;
    mTotalHeatFlux      = 0.0;

    // Compute heat fluxes with neighbor particles
    for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
      if (mNeighbourElements[i] == NULL) continue;
      mNeighbor_p   = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(mNeighbourElements[i]);
      mNeighborType = PARTICLE_NEIGHBOR;
      ComputeHeatFluxWithNeighbor(r_process_info);
    }

    // Compute heat fluxes with neighbor walls
    for (unsigned int i = 0; i < mNeighbourRigidFaces.size(); i++) {
      if (mNeighbourRigidFaces[i] == NULL) continue;
      mNeighbor_w   = dynamic_cast<DEMWall*>(mNeighbourRigidFaces[i]);
      mNeighborType = WALL_NEIGHBOR;
      ComputeHeatFluxWithNeighbor(r_process_info);
    }

    // Compute convection with surrounding fluid
    if (this->Is(DEMFlags::HAS_CONVECTION))
      ComputeConvectiveHeatFlux(r_process_info);

    // Sum up contributions
    mTotalHeatFlux = mConductiveHeatFlux + mConvectiveHeatFlux + mRadiativeHeatFlux + mPrescribedHeatFlux;

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Finalization methods

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::FinalizeSolutionStep(const ProcessInfo& r_process_info) {
    TBaseElement::FinalizeSolutionStep(r_process_info);
    UpdateTemperature(r_process_info);
    mPreviousTemperature = GetParticleTemperature();
    SetParticleHeatFlux(mTotalHeatFlux);
  }

  //=====================================================================================================================================================================================
  // Update methods

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::UpdateTemperature(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    if (!this->Is(DEMFlags::HAS_FIXED_TEMPERATURE) && !this->Is(DEMFlags::IS_ADIABATIC)) {
      // Compute new temperature
      double thermal_inertia = GetMass() * GetParticleHeatCapacity();
      double temp_increment  = mTotalHeatFlux / thermal_inertia * r_process_info[DELTA_TIME];
      double temp_new        = GetParticleTemperature() + temp_increment;
      
      // Set new temperature
      SetParticleTemperature(temp_new);
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::UpdateTemperatureDependentRadius(const ProcessInfo& r_process_info) {
    //double this_temp      = GetParticleTemperature();
    //double fluid_temp     = r_process_info[FLUID_TEMPERATURE];
    //double relative_temp  = this_temp - fluid_temp; // temp in Kelvin
    //double thermal_alpha  = GetProperties()[THERMAL_EXPANSION_COEFFICIENT];
    //double updated_radius = GetRadius() * (1 + thermal_alpha * relative_temp);
    //SetRadius(updated_radius);
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion(const ProcessInfo& r_process_info,
                                                                                                              double& thermalDeltDisp,
                                                                                                              double& thermalRelVel,
                                                                                                              ThermalSphericParticle<TBaseElement>* element2) {
    //thermalRelVel                      = 0;
    //double temperature                 = GetParticleTemperature();
    //double other_temperature           = element2->GetParticleTemperature();
    //double previous_temperature        = mPreviousTemperature;
    //double previous_other_temperature  = element2->mPreviousTemperature;
    //double thermal_alpha               = GetProperties()[THERMAL_EXPANSION_COEFFICIENT];
    //double updated_radius              = GetRadius();
    //double updated_other_radius        = element2->GetRadius();
    //double dt                          = r_process_info[DELTA_TIME];
    //double temperature_increment_elem1 = temperature - previous_temperature;
    //double temperature_increment_elem2 = other_temperature - previous_other_temperature;
    //thermalDeltDisp                    = updated_radius * thermal_alpha * temperature_increment_elem1;
    //thermalDeltDisp                    = thermalDeltDisp + updated_other_radius * thermal_alpha * temperature_increment_elem2;
    //thermalRelVel                      = updated_radius * thermal_alpha * temperature_increment_elem1 / dt;
    //thermalRelVel                      = thermalRelVel + updated_other_radius * thermal_alpha * temperature_increment_elem2 / dt;
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(const ProcessInfo& r_process_info,
                                                                                                            double DeltDisp[3], //IN GLOBAL AXES
                                                                                                            double RelVel[3],   //IN GLOBAL AXES
                                                                                                            double OldLocalCoordSystem[3][3],
                                                                                                            double LocalCoordSystem[3][3],
                                                                                                            SphericParticle* neighbor_iterator) {
    //double thermalDeltDisp = 0;
    //double thermalRelVel   = 0;
    //ThermalSphericParticle<TBaseElement>* thermal_neighbor_iterator = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(neighbor_iterator);
    //UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion(r_process_info, thermalDeltDisp, thermalRelVel, thermal_neighbor_iterator);
    //double LocalRelVel[3] = {0.0};
    //GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel); //TODO: can we do this in global axes directly?
    ////LocalRelVel[2] -= thermalRelVel;
    //GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalRelVel, RelVel);
  }

  //=====================================================================================================================================================================================
  // Heat fluxes computation

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeHeatFluxWithNeighbor(const ProcessInfo& r_process_info) {
    KRATOS_TRY

  // Check if neighbor is adiabatic
  if (CheckAdiabaticNeighbor())
    return;

  // Set distance between particles centroids
  SetDistanceToNeighbor();

  // Heat transfer mechanisms
  if (this->Is(DEMFlags::HAS_DIRECT_CONDUCTION))   ComputeDirectConductionHeatFlux(r_process_info);
  if (this->Is(DEMFlags::HAS_INDIRECT_CONDUCTION)) ComputeIndirectConductionHeatFlux(r_process_info);
  if (this->Is(DEMFlags::HAS_RADIATION))           ComputeRadiativeHeatFlux(r_process_info);

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeDirectConductionHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Check for contact
    if (!CheckHeatTransferDistance(0.0))
      return;

    // Compute heat flux according to selected model
    std::string model = r_process_info[DIRECT_CONDUCTION_MODEL];

    if      (model.compare("batchelor_obrien") == 0) mConductiveHeatFlux += DirectConductionBatchelorOBrien(r_process_info);
    else if (model.compare("thermal_pipe")     == 0) mConductiveHeatFlux += DirectConductionThermalPipe(r_process_info);
    else if (model.compare("collisional")      == 0) mConductiveHeatFlux += DirectConductionCollisional(r_process_info);

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeIndirectConductionHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Compute heat flux according to selected model
    std::string model = r_process_info[INDIRECT_CONDUCTION_MODEL];

    if      (model.compare("surrounding_layer") == 0) mConductiveHeatFlux += IndirectConductionSurroundingLayer(r_process_info);
    else if (model.compare("voronoi_a")         == 0) mConductiveHeatFlux += IndirectConductionVoronoiA(r_process_info);
    else if (model.compare("voronoi_b")         == 0) mConductiveHeatFlux += IndirectConductionVoronoiB(r_process_info);

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeRadiativeHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // TODO: radiation with walls not yet implemented
    if (mNeighborType == WALL_NEIGHBOR)
      return;

    // Compute heat flux according to selected model
    std::string model = r_process_info[RADIATION_MODEL];

    if      (model.compare("continuum_zhou")   == 0) mRadiativeHeatFlux += RadiationContinuumZhou(r_process_info);
    else if (model.compare("continuum_krause") == 0) mRadiativeHeatFlux += RadiationContinuumKrause(r_process_info);

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeConvectiveHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY
    
    double surface_area       = GetParticleSurfaceArea();
    double char_length        = GetParticleCharacteristicLength();
    double fluid_conductivity = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    double temp_grad          = r_process_info[FLUID_TEMPERATURE] - GetParticleTemperature();

    // Compute Nusselt number according to selected model
    double Nu = 0.0;
    std::string model = r_process_info[CONVECTION_MODEL];
    if      (model.compare("sphere_hanz_marshall") == 0) Nu = ConvectionHanzMarshall(r_process_info);
    else if (model.compare("sphere_whitaker")      == 0) Nu = ConvectionWhitaker(r_process_info);

    // Compute heat flux
    mConvectiveHeatFlux += (Nu * fluid_conductivity / char_length) * surface_area * temp_grad;

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Heat transfer models

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::DirectConductionBatchelorOBrien(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    double keff      = ComputeEffectiveConductivity();
    double Rc        = ComputeContactRadius();
    double temp_grad = GetNeighborTemperature() - GetParticleTemperature();

    return 4.0 * keff * Rc * temp_grad;

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::DirectConductionThermalPipe(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    double kavg      = ComputeAverageConductivity();
    double Rc        = ComputeContactRadius();
    double temp_grad = GetNeighborTemperature() - GetParticleTemperature();

    return kavg * (Globals::Pi * Rc * Rc) * temp_grad / mNeighborDistance;

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::DirectConductionCollisional(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Check if collision time is smaller than expected value, otherwise use static model (batchelor_obrien)
    // TODO: track collision time and save impact velocity
    double col_time = 0.0;
    double impact_normal_velocity = 0.0;
    double col_time_max = ComputeMaxCollisionTime();

    if (col_time < col_time_max && impact_normal_velocity != 0.0) {
      double temp_grad = GetNeighborTemperature() - GetParticleTemperature();
      double Rc_max    = ComputeMaxContactRadius();
      double Fo        = ComputeFourierNumber();

      double a1 = GetDensity() * GetParticleHeatCapacity();
      double a2 = GetNeighborDensity() * GetNeighborHeatCapacity();
      double b1 = a1 * GetParticleConductivity();
      double b2 = a2 * GetNeighborConductivity();
      double c  = a1 / a2;

      double C1 = -2.300 * c * c +  8.909 * c - 4.235;
      double C2 =  8.169 * c * c - 33.770 * c + 24.885;
      double C3 = -5.758 * c * c + 24.464 * c - 20.511;

      double C_coeff = 0.435 * (sqrt(C2 * C2 - 4 * C1 * (C3 - Fo)) - C2) / C1;

      return C_coeff * Globals::Pi * Rc_max * Rc_max * pow(col_time_max,-1/2) * temp_grad / (pow(b1,-1/2) + pow(b2,-1/2));
    }
    else {
      return DirectConductionBatchelorOBrien(r_process_info);
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::IndirectConductionSurroundingLayer(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Check if particles are close enough
    double layer = r_process_info[FLUID_LAYER_THICKNESS];
    if (!CheckHeatTransferDistance(layer))
      return 0.0;

    // Compute heat transfer coefficient
    // Assumption: neighbor wall is treated as a particle with the same radius
    double particle_radius    = GetRadius();
    double neighbor_radius    = (mNeighborType == PARTICLE_NEIGHBOR) ? mNeighbor_p->GetRadius() : 0.0;
    double fluid_conductivity = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    double h = 0.0;

    if (particle_radius == neighbor_radius || mNeighborType == WALL_NEIGHBOR) {
      double min_dist = r_process_info[MIN_CONDUCTION_DISTANCE];
      double a = (mNeighborDistance - 2.0 * particle_radius) / particle_radius;
      double r_in, r_out, b, c;

      if (mNeighborDistance > 2.0 * particle_radius + min_dist)
        r_in = 0.0;
      else
        r_in = sqrt(1.0 - pow(min_dist / particle_radius - a - 1.0, 2));

      if (a > sqrt(pow((particle_radius + (layer * particle_radius)) / particle_radius, 2) - 1.0) - 1.0)
        r_out = sqrt(pow((particle_radius + (layer * particle_radius)) / particle_radius, 2) - (a + 1.0) * (a + 1.0));
      else
        r_out = 1.0;

      b = sqrt(1.0 - r_out * r_out);
      c = sqrt(1.0 - r_in  * r_in);

      // Heat transfer coefficient from analytica solution of the integral expression
      h = 2.0 * Globals::Pi * fluid_conductivity * particle_radius * ((a + 1.0) * log(abs((b - a - 1.0) / (a - c + 1.0))) + b - c);
    }
    else {
      // Compute lower limit of integral (contact radius)
      double low_lim = ComputeContactRadius();

      // Compute upper limit of integral
      double r_min = std::min(particle_radius, neighbor_radius);
      double r_max = std::max(particle_radius, neighbor_radius);
      double param = pow((r_max + (layer * r_max)), 2);
      double upp_lim;

      if (mNeighborDistance <= sqrt(param - r_min * r_min))
        upp_lim = r_min;
      else
        upp_lim = sqrt(param - pow(((param - r_min * r_min + mNeighborDistance * mNeighborDistance) / (2.0 * mNeighborDistance)), 2));

      // Build struct of integration parameters
      struct IntegrandParams params;
      params.r_process_info = r_process_info;
      params.r1 = particle_radius;
      params.r2 = neighbor_radius;

      // Heat transfer coefficient from integral expression solved numerically
      h = fluid_conductivity * AdaptiveSimpsonIntegration(low_lim, upp_lim, &EvalIntegrandSurrLayer(params));
    }

    // Compute heat flux
    return h * GetNeighborTemperature() - GetParticleTemperature();

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::IndirectConductionVoronoiA(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Check if particles are close enough
    if (!CheckHeatTransferDistance(r_process_info[MAX_CONDUCTION_DISTANCE]))
      return 0.0;

    // Compute lower limit of integral (contact radius)
    double particle_radius = GetRadius();
    double neighbor_radius = (mNeighborType == PARTICLE_NEIGHBOR) ? mNeighbor_p->GetRadius() : 0.0;
    double low_lim = ComputeContactRadius();

    // Compute voronoi edge radius from porosity
    // TODO: currently, global average porosity is an input
    double porosity = r_process_info[PRESCRIBED_GLOBAL_POROSITY];
    double rij = 0.56 * particle_radius / pow((1.0 - porosity), 1/3);

    // Compute heat transfer coefficient
    // Assumption: neighbor wall is treated as a particle with the same radius
    double h = 0.0;

    if (particle_radius == neighbor_radius || mNeighborType == WALL_NEIGHBOR) {
      double keff = ComputeEffectiveConductivity();
      double upp_lim = particle_radius * rij / sqrt(rij * rij + mNeighborDistance * mNeighborDistance / 4.0);

      // Build struct of integration parameters
      struct IntegrandParams params;
      params.r_process_info = r_process_info;
      params.keff = ComputeEffectiveConductivity();
      params.r1   = particle_radius;
      params.rij  = rij;

      // Heat transfer coefficient from integral expression solved numerically
      h = AdaptiveSimpsonIntegration(low_lim, upp_lim, &EvalIntegrandVoronoiMono(params));
    }
    else {
      double D1, D2, rij_, upp_lim;

      if (mNeighborDistance < particle_radius + neighbor_radius)
        D1 = sqrt(particle_radius * particle_radius - low_lim * low_lim);
      else
        D1 = (particle_radius * particle_radius - neighbor_radius * neighbor_radius + mNeighborDistance * mNeighborDistance) / (2 * mNeighborDistance);

      D2 = mNeighborDistance - D1;

      if (particle_radius <= neighbor_radius)
        upp_lim = particle_radius * rij / sqrt(rij * rij + D1 * D1);
      else
        upp_lim = neighbor_radius * rij / sqrt(rij * rij + D2 * D2);

      rij_ = D2 * upp_lim / sqrt(neighbor_radius * neighbor_radius - upp_lim * upp_lim);

      // Build struct of integration parameters
      struct IntegrandParams params;
      params.r_process_info = r_process_info;
      params.r1   = particle_radius;
      params.r2   = neighbor_radius;
      params.rij  = rij;
      params.rij_ = rij_;
      params.D1   = D1;
      params.D2   = D2;
      params.k1   = GetParticleConductivity();
      params.k2   = GetNeighborConductivity();

      // Heat transfer coefficient from integral expression solved numerically
      h = AdaptiveSimpsonIntegration(low_lim, upp_lim, &EvalIntegrandVoronoiMulti(params));
    }

    // Compute heat flux
    return h * GetNeighborTemperature() - GetParticleTemperature();

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::IndirectConductionVoronoiB(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Check if particles are close enough
    if (!CheckHeatTransferDistance(r_process_info[MAX_CONDUCTION_DISTANCE]))
      return 0.0;

    // Get parameters
    double Rc                    = ComputeContactRadius();
    double particle_radius       = GetRadius();
    double neighbor_radius       = (mNeighborType == PARTICLE_NEIGHBOR) ? mNeighbor_p->GetRadius() : 0.0;
    double particle_conductivity = GetParticleConductivity();
    double neighbor_conductivity = GetNeighborConductivity();
    double fluid_conductivity    = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    double core                  = r_process_info[ISOTHERMAL_CORE_RADIUS];

    // Compute voronoi edge radius from porosity
    // TODO: currently, global average porosity is an input
    double porosity = r_process_info[PRESCRIBED_GLOBAL_POROSITY];
    double rij = 0.56 * particle_radius / pow((1.0 - porosity), 1/3);

    // Compute heat transfer coefficient
    // Assumption: neighbor wall is treated as a particle with the same radius
    double h = 0.0;

    if (particle_radius == neighbor_radius || mNeighborType == WALL_NEIGHBOR) {
      double keff = ComputeEffectiveConductivity();
      double D    = mNeighborDistance / 2.0;
      double a    = (1.0 / core - 1.0 / particle_radius) / (2.0 * keff) + 1.0 / (fluid_conductivity * particle_radius);
      double b    = 1.0 / (fluid_conductivity * D);
      double c0   = D / sqrt(rij * rij + D * D);
      double c1   = D / sqrt(Rc * Rc + D * D);
      double f    = (a - b * c0) / (a - b * c1);
      double ln   = 0.0;
      if (f > 0.0)
        ln = log(f);

      // Heat transfer coefficient
      h = Globals::Pi * ln / b;
    }
    else {
      double An = Globals::Pi * rij * rij;  // Compute area of neighboring voronoi cells

      double gamma1 = particle_radius / mNeighborDistance;
      double gamma2 = neighbor_radius / mNeighborDistance;
      double dgamma = gamma2 - gamma1;

      double A = (particle_conductivity + fluid_conductivity * (1.0 / core - 1.0)) / (particle_conductivity * gamma1);
      double B = (neighbor_conductivity + fluid_conductivity * (1.0 / core - 1.0)) / (neighbor_conductivity * gamma2);

      double lambda = (1.0 + dgamma * A) * (1.0 - dgamma * B);

      double delmax = 0.5 * (sqrt((4.0 * An) / (Globals::Pi * mNeighborDistance * mNeighborDistance * (1.0 - dgamma * dgamma)) + 1.0) - dgamma);
      double delmin = 0.5 * (sqrt((4.0 * Rc * Rc) / (mNeighborDistance * mNeighborDistance * (1.0 - dgamma * dgamma)) + 1.0) - dgamma);

      double Xmax = ((A + B) * delmax + dgamma * B - 1.0) / sqrt(fabs(lambda));
      double Xmin = ((A + B) * delmin + dgamma * B - 1.0) / sqrt(fabs(lambda));

      double Y1 = (Xmax - Xmin) / (1.0 - Xmax * Xmin);
      double Y2 = (Xmax - Xmin) / (1.0 + Xmax * Xmin);

      // Heat transfer coefficient
      if (lambda > 0.0)
        h = Globals::Pi * fluid_conductivity * mNeighborDistance * (1.0 - dgamma * dgamma) * log(fabs((1.0 - Y1) / (1.0 + Y1))) / (2.0 * sqrt(fabs(lambda)));
      else if (lambda < 0.0)
        h = Globals::Pi * fluid_conductivity * mNeighborDistance * (1.0 - dgamma * dgamma) * atan(Y2) / (2.0 * sqrt(fabs(lambda)));
      else
        h = Globals::Pi * fluid_conductivity * mNeighborDistance * (1.0 - dgamma * dgamma) * (1.0 / delmin - 1.0 / delmax) / (A + B);
    }

    // Compute heat flux
    return h * GetNeighborTemperature() - GetParticleTemperature();

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::RadiationContinuumZhou(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // TODO: Needs the porosity
    return 0.0;

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::RadiationContinuumKrause(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // TODO: Needs to be computed outsie neighbors loop
    return 0.0;

    //// Get particle properties
    //double this_radius = GetRadius();
    //double this_area = 4 * Globals::Pi * this_radius * this_radius;
    //double this_emissivity = mThermalEmissivity;
    //double this_temp = GetParticleTemperature();

    //// Initialize parameters
    //double num = 0.0;
    //double den = 0.0;

    //// Loop over neighbor particles
    //// (Assumption: walls not considered)
    //for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
    //  if (mNeighbourElements[i] == NULL) continue;
    //  ThermalSphericParticle<TBaseElement>* neighbor_iterator = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(mNeighbourElements[i]);

    //  // Check if neighbor is adiabatic
    //  if (neighbor_iterator->Is(DEMFlags::IS_ADIABATIC))
    //    continue;

    //  // Compute direction and distance between centroids
    //  array_1d<double, 3> direction;
    //  direction[0] = GetGeometry()[0].Coordinates()[0] - neighbor_iterator->GetGeometry()[0].Coordinates()[0];
    //  direction[1] = GetGeometry()[0].Coordinates()[1] - neighbor_iterator->GetGeometry()[0].Coordinates()[1];
    //  direction[2] = GetGeometry()[0].Coordinates()[2] - neighbor_iterator->GetGeometry()[0].Coordinates()[2];

    //  double distance = DEM_MODULUS_3(direction);

    //  // Check if particles are close enough
    //  // (Assumption: radiation influence factor applied to the maximum radius)
    //  double other_radius = neighbor_iterator->GetRadius();
    //  if (distance > r_process_info[RADIATION_RADIUS] * std::max(this_radius, other_radius))
    //    continue;

    //  // Get particle properties
    //  double other_area = 4 * Globals::Pi * other_radius * other_radius;
    //  double other_emissivity = neighbor_iterator->mThermalEmissivity;
    //  double other_temp = neighbor_iterator->GetParticleTemperature();

    //  // Update parameters
    //  den += STEFAN_BOLTZMANN * other_emissivity * other_area / 2;
    //  num += den * pow(other_temp, 4);
    //}

    //// Averaged environment temperature
    //double env_temp = pow(num / den, 1 / 4);

    //// Compute heat flux
    //mRadiativeHeatFlux += this_emissivity * STEFAN_BOLTZMANN * this_area * (pow(env_temp, 4) - pow(this_temp, 4));

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::ConvectionHanzMarshall(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    double Pr = ComputePrandtlNumber(r_process_info);
    double Re = ComputeReynoldNumber(r_process_info);

    return 2.0 + 0.6 * pow(Re,1/2) * pow(Pr,1/3);

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::ConvectionWhitaker(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    double Pr = ComputePrandtlNumber(r_process_info);
    double Re = ComputeReynoldNumber(r_process_info);

    // Assumption: temperature-dependent viscosity at particle surface is negleted
    return 2.0 + (0.4 * pow(Re,1/2) + 0.06 * pow(Re,2/3)) * pow(Pr,2/5);

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Auxiliary computations

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeAddedSearchDistance(const ProcessInfo& r_process_info, double& added_search_distance) {
    KRATOS_TRY

    if (this->Is(DEMFlags::HAS_INDIRECT_CONDUCTION)){
      std::string model = r_process_info[INDIRECT_CONDUCTION_MODEL];
      if (model.compare("surrounding_layer") == 0) {
        double model_search_distance = GetRadius() * r_process_info[FLUID_LAYER_THICKNESS];
        added_search_distance = std::max(added_search_distance, model_search_distance);
      }
      else if (model.compare("voronoi_a") == 0 ||
               model.compare("voronoi_b") == 0) {
        double model_search_distance = GetRadius() * r_process_info[MAX_CONDUCTION_DISTANCE];
        added_search_distance = std::max(added_search_distance, model_search_distance);
      }
    }

    if (this->Is(DEMFlags::HAS_RADIATION)) {
      std::string model = r_process_info[RADIATION_MODEL];
      if (model.compare("sphere_hanz_marshall") == 0 ||
          model.compare("sphere_whitaker")      == 0) {
        double model_search_distance = GetRadius() * (r_process_info[RADIATION_RADIUS] - 1);
        added_search_distance = std::max(added_search_distance, model_search_distance);
      }
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::ComputePrandtlNumber(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    double fluid_viscosity    = r_process_info[FLUID_VISCOSITY];
    double fluid_heatcapacity = r_process_info[FLUID_HEAT_CAPACITY];
    double fluid_conductivity = r_process_info[FLUID_THERMAL_CONDUCTIVITY];

    return fluid_viscosity * fluid_heatcapacity / fluid_conductivity;

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::ComputeReynoldNumber(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    double char_length     = GetParticleCharacteristicLength();
    double rel_velocity    = ComputeFluidRelativeVelocity(r_process_info);
    double fluid_density   = r_process_info[FLUID_DENSITY];
    double fluid_viscosity = r_process_info[FLUID_VISCOSITY];

    return fluid_density * char_length * rel_velocity / fluid_viscosity;

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::ComputeFluidRelativeVelocity(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    array_1d<double, 3> particle_velocity = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    array_1d<double, 3> rel_velocity      = r_process_info[FLUID_VELOCITY];

    for (unsigned int i = 0; i < rel_velocity.size(); i++)
      rel_velocity[i] -= particle_velocity[i];

    return DEM_MODULUS_3(rel_velocity);

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::AdaptiveSimpsonIntegration(double a, double b, double (*evalIntegrand)(IntegrandParams params)) {
    KRATOS_TRY

    // Initialization
    params.position = a;
    double fa = evalIntegrand(params);

    params.position = b;
    double fb = evalIntegrand(params);

    params.position = (a + b) / 2.0;
    double fc = evalIntegrand(params);

    // Get tolerance
    double tol = params.r_process_info[INTEGRAL_TOLERANCE];
    constexpr double eps = std::numeric_limits<double>::epsilon();
    if (tol < 10.0 * eps) tol = 10.0 * eps;

    // Solve integral recursively with adaptive Simpson quadrature
    return RecursiveSimpsonIntegration(a, b, fa, fb, fc, tol, evalIntegrand);

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::RecursiveSimpsonIntegration(double a, double b, double fa, double fb, double fc, double tol, double (*evalIntegrand)(IntegrandParams params)) {
    KRATOS_TRY

    // TODO: in order to catch possible erros that can occur in singularities,
    //       add a min value for subdivision size (to contain machine representable point) and a max number of function evaluation (+- 10000).

    double c  = (a + b) / 2.0;

    params.position = (a + c) / 2.0;
    double fd = evalIntegrand(params);

    params.position = (c + b) / 2.0;
    double fe = evalIntegrand(params);

    double I1 = (b - a) / 6.0  * (fa + 4.0 * fc + fb);
    double I2 = (b - a) / 12.0 * (fa + 4.0 * fd + 2.0 * fc + 4.0 * fe + fb);

    if (fabs(I2-I1) <= tol) {
      return I2 + (I2 - I1) / 15.0;
    }
    else { // sub-divide interval recursively
      double Ia = RecursiveSimpsonIntegration(a, c, fa, fc, fd, tol, evalIntegrand);
      double Ib = RecursiveSimpsonIntegration(c, b, fc, fb, fe, tol, evalIntegrand);
      return Ia + Ib;
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::EvalIntegrandSurrLayer(IntegrandParams params) {
    KRATOS_TRY

    double r  = params.position;
    double r1 = params.r1;
    double r2 = params.r2;
    double d  = params.r_process_info[MIN_CONDUCTION_DISTANCE];

    return 2.0 * Globals::Pi * r / std::max(d, mNeighborDistance - sqrt(r1 * r1 - r * r) - sqrt(r2 * r2 - r * r));

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::EvalIntegrandVoronoiMono(IntegrandParams params) {
    KRATOS_TRY

    double r    = params.position;
    double rp   = params.r1;
    double rij  = params.rij;
    double keff = params.keff;
    double kf   = params.r_process_info[FLUID_THERMAL_CONDUCTIVITY];

    return 2.0 * Globals::Pi * r / ((sqrt(rp * rp - r * r) - r * mNeighborDistance / (2.0 * rij)) / keff + 2.0 * (mNeighborDistance / 2.0 - sqrt(rp * rp - r * r)) / kf);

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::EvalIntegrandVoronoiMulti(IntegrandParams params) {
    KRATOS_TRY

    double r     = params.position;
    double r1    = params.r1;
    double r2    = params.r2;
    double rij   = params.rij;
    double rij_  = params.rij_;
    double D1    = params.D1;
    double D2    = params.D2;
    double k1    = params.k1;
    double k2    = params.k2;
    double kf    = params.r_process_info[FLUID_THERMAL_CONDUCTIVITY];

    return 2.0 * Globals::Pi * r / ((sqrt(r1 * r1 - r * r) - r * D1 / rij) / k1 + (sqrt(r2 * r2 - r * r) - r * D2 / rij_) / k2 + (mNeighborDistance - sqrt(r1 * r1 - r * r) - sqrt(r2 * r2 - r * r)) / kf);

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Neighbor interaction computations

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeContactArea(const double rmin, double indentation, double& calculation_area) {
    calculation_area = Globals::Pi*rmin*rmin;
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::SetDistanceToNeighbor() {
    KRATOS_TRY

    if (mNeighborType == PARTICLE_NEIGHBOR) {
      array_1d<double, 3> direction;
      direction[0] = GetGeometry()[0].Coordinates()[0] - mNeighbor_p->GetGeometry()[0].Coordinates()[0];
      direction[1] = GetGeometry()[0].Coordinates()[1] - mNeighbor_p->GetGeometry()[0].Coordinates()[1];
      direction[2] = GetGeometry()[0].Coordinates()[2] - mNeighbor_p->GetGeometry()[0].Coordinates()[2];
      mNeighborDistance = DEM_MODULUS_3(direction);
    }
    else if (mNeighborType == WALL_NEIGHBOR) {
      // Stolen from ComputeBallToRigidFaceContactForce in spheric_particle
      double dummy1[3][3];
      DEM_SET_COMPONENTS_TO_ZERO_3x3(dummy1);
      array_1d<double, 4> dummy2 = this->mContactConditionWeights[i];
      array_1d<double, 3> dummy3 = ZeroVector(3);
      array_1d<double, 3> dummy4 = ZeroVector(3);
      int dummy5 = -1;
      double distance = 0.0;
      mNeighbor_w->ComputeConditionRelativeData(i, this, dummy1, distance, dummy2, dummy3, dummy4, dummy5);
      mNeighborDistance = distance;

      // Stolen from ComputeBallToRigidFaceContactForce in spheric_particle
      //array_1d<double, 3> cond_to_me_vect;
      //array_1d<double, 3> wall_coordinates = mNeighbor_w->GetGeometry().Center();
      //noalias(cond_to_me_vect) = GetGeometry()[0].Coordinates() - wall_coordinates;
      //mNeighborDistance = DEM_MODULUS_3(cond_to_me_vect);
    }
    else {
      mNeighborDistance = 0.0;
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  bool ThermalSphericParticle<TBaseElement>::CheckHeatTransferDistance(const double radius_factor) {
    KRATOS_TRY

    double particle_radius = GetRadius();
    double neighbor_radius = (mNeighborType == PARTICLE_NEIGHBOR) ? mNeighbor_p->GetRadius() : 0.0;
    double separation      = mNeighborDistance - particle_radius - neighbor_radius

    // Assumption: radius_factor applies to the larger radius
    return (separation < radius_factor * std::max(particle_radius, neighbor_radius));

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  bool ThermalSphericParticle<TBaseElement>::CheckAdiabaticNeighbor(void) {
    KRATOS_TRY

    return ((mNeighborType == PARTICLE_NEIGHBOR && mNeighbor_p->Is(DEMFlags::IS_ADIABATIC)) ||
            (mNeighborType == WALL_NEIGHBOR     && mNeighbor_w->Is(DEMFlags::IS_ADIABATIC)));

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::ComputeFourierNumber() {
    KRATOS_TRY

    double col_time_max = ComputeMaxCollisionTime();
    double Rc_max       = ComputeMaxContactRadius();

    double Fo_particle = GetParticleConductivity() * col_time_max / (GetDensity() * GetParticleHeatCapacity() * Rc_max * Rc_max);
    double Fo_neighbor;

    if (mNeighborType == PARTICLE_NEIGHBOR)
      Fo_neighbor = GetNeighborConductivity() * col_time_max / (GetNeighborDensity() * GetNeighborHeatCapacity() * Rc_max * Rc_max);
    else if (mNeighborType == WALL_NEIGHBOR)
      Fo_neighbor = 0.0;

    // Assumption: average of both particles
    return (Fo_particle + Fo_neighbor) / 2.0;

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::ComputeMaxCollisionTime() {
    KRATOS_TRY

    // TODO: save impact normal velocity
    double impact_normal_velocity = 0.0;
    double eff_radius = ComputeEffectiveRadius();
    double eff_mass   = ComputeEffectiveMass();
    double eff_young  = ComputeEffectiveYoung();

    return 2.87 * pow(eff_mass * eff_mass / (eff_radius * eff_young * eff_young * impact_normal_velocity), 1/5);

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::ComputeMaxContactRadius() {
    KRATOS_TRY

    // TODO: save impact normal velocity
    double impact_normal_velocity = 0.0;
    double eff_radius = ComputeEffectiveRadius();
    double eff_mass   = ComputeEffectiveMass();
    double eff_young  = ComputeEffectiveYoung();

    return pow(15.0 * eff_radius * eff_mass * impact_normal_velocity * impact_normal_velocity / (16.0 * eff_young), 1/5);

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::ComputeContactRadius() {
    KRATOS_TRY
    
    double Rc = 0.0;

    if (mNeighborType == PARTICLE_NEIGHBOR) {
      double r1 = GetRadius();
      double r2 = mNeighbor_p->GetRadius();
      if (mNeighborDistance < r1 + r2)
        Rc = sqrt(fabs(r1 * r1 - pow(((r1 * r1 - r2 * r2 + mNeighborDistance * mNeighborDistance) / (2.0 * mNeighborDistance)), 2)));
      else
        Rc = 0.0;
    }
    else if (mNeighborType == WALL_NEIGHBOR) {
      double r = GetRadius();
      double ident = r - mNeighborDistance;
      if (ident > 0.0)
        Rc = sqrt(ident * (2.0 * r - ident));
      else
        Rc = 0.0;
    }

    return Rc;

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::ComputeEffectiveRadius() {
    KRATOS_TRY

    if (mNeighborType == PARTICLE_NEIGHBOR) {
      double particle_radius = GetRadius();
      double neighbor_radius = mNeighbor_p->GetRadius();
      return particle_radius * neighbor_radius / (particle_radius + neighbor_radius);
    }
    else if (mNeighborType == WALL_NEIGHBOR) {
      return GetRadius();
    }
    else {
      return 0.0;
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::ComputeEffectiveMass() {
    KRATOS_TRY

    if (mNeighborType == PARTICLE_NEIGHBOR) {
      double particle_mass = GetMass();
      double neighbor_mass = mNeighbor_p->GetMass();
      return particle_mass * neighbor_mass / (particle_mass + neighbor_mass);
    }
    else if (mNeighborType == WALL_NEIGHBOR) {
      return GetMass();
    }
    else {
      return 0.0;
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::ComputeEffectiveYoung() {
    KRATOS_TRY

    double particle_young   = GetYoung();
    double particle_poisson = GetPoisson();
    double neighbor_young   = GetNeighborYoung();
    double neighbor_poisson = GetNeighborPoisson();

    return 1 / ((1 - particle_poisson * particle_poisson) / particle_young + (1 - neighbor_poisson * neighbor_poisson) / neighbor_young);

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::ComputeEffectiveConductivity() {
    KRATOS_TRY

    if (mNeighborType == PARTICLE_NEIGHBOR) {
      double particle_conductivity = GetParticleConductivity();
      double neighbor_conductivity = GetNeighborConductivity();
      return particle_conductivity * neighbor_conductivity / (particle_conductivity + neighbor_conductivity);
    }
    else if (mNeighborType == WALL_NEIGHBOR) {
      return GetParticleConductivity();
    }
    else {
      return 0.0;
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::ComputeAverageConductivity() {
    KRATOS_TRY

    if (mNeighborType == PARTICLE_NEIGHBOR) {
      double r1 = GetRadius();
      double r2 = mNeighbor_p->GetRadius();
      double particle_conductivity = GetParticleConductivity();
      double neighbor_conductivity = GetNeighborConductivity();
      return (r1 + r2) / (r1 / particle_conductivity + r2 / neighbor_conductivity);
    }
    else if (mNeighborType == WALL_NEIGHBOR) {
      // Assumption: average conductivity considers particle only
      return GetParticleConductivity();
    }
    else {
      return 0.0;
    }

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Get/Set methods

  template <class TBaseElement>
  const double& ThermalSphericParticle<TBaseElement>::GetParticleTemperature() {
    return GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::GetParticleSurfaceArea() {
    return 4 * Globals::Pi * GetRadius() * GetRadius();
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::GetParticleCharacteristicLength() {
    return 2 * GetRadius();
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::GetParticleConductivity() {
    return GetProperties()[THERMAL_CONDUCTIVITY];
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::GetParticleHeatCapacity() {
    return GetProperties()[SPECIFIC_HEAT];
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::GetParticleEmissivity() {
    return GetProperties()[EMISSIVITY];
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::GetWallTemperature() {
    // Assumption: wall temperature is the average of its nodes
    double wall_temp = 0.0;
    double n_nodes = mNeighbor_w->GetGeometry().size();
    for (unsigned int i = 0; i < n_nodes; i++)
      wall_temp += mNeighbor_w->GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE);
    return wall_temp / n_nodes;
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::GetNeighborTemperature() {
    if (mNeighborType == PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleTemperature();
    else if (mNeighborType == WALL_NEIGHBOR)
      return mNeighbor_w->GetWallTemperature();
    else
      return 0.0;
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::GetNeighborDensity() {
    if (mNeighborType == PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetDensity();
    else if (mNeighborType == WALL_NEIGHBOR)
      return mNeighbor_w->GetProperties()[DENSITY];
    else
      return 0.0;
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::GetNeighborYoung() {
    if (mNeighborType == PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetYoung();
    else if (mNeighborType == WALL_NEIGHBOR)
      return mNeighbor_w->GetProperties()[YOUNG_MODULUS];
    else
      return 0.0;
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::GetNeighborPoisson() {
    if (mNeighborType == PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetPoisson();
    else if (mNeighborType == WALL_NEIGHBOR)
      return mNeighbor_w->GetProperties()[POISSON_RATIO];
    else
      return 0.0;
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::GetNeighborConductivity() {
    if (mNeighborType == PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleConductivity();
    else if (mNeighborType == WALL_NEIGHBOR)
      return mNeighbor_w->GetProperties()[THERMAL_CONDUCTIVITY];
    else
      return 0.0;
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::GetNeighborHeatCapacity() {
    if (mNeighborType == PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleHeatCapacity();
    else if (mNeighborType == WALL_NEIGHBOR)
      return mNeighbor_w->GetProperties()[SPECIFIC_HEAT];
    else
      return 0.0;
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::SetParticleTemperature(const double temperature) {
    GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE) = temperature;
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::SetParticleHeatFlux(const double heat_flux) {
    GetGeometry()[0].GetSolutionStepValue(HEATFLUX) = heat_flux;
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::SetParticlePrescribedHeatFlux(const double heat_flux) {
    mPrescribedHeatFlux = heat_flux;
  }

  // Explicit Instantiation
  template class ThermalSphericParticle<SphericParticle>;
  template class ThermalSphericParticle<SphericContinuumParticle>;

} // namespace Kratos
