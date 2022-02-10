//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics ThermalDEM Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Rafael Rangel (rrangel@cimne.upc.edu)
//

/*
* Important tags:
*   attention;
*   assumption;
*   todo;
*/

// System includes

// External includes

// Project includes
#include "thermal_spheric_particle.h"

namespace Kratos
{
  //=====================================================================================================================================================================================
  // Constructor/Destructor methods

  ThermalSphericParticle::ThermalSphericParticle():SphericParticle() {
    mpThermalIntegrationScheme = NULL;
  }

  ThermalSphericParticle::ThermalSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry):SphericParticle(NewId, pGeometry) {
    mpThermalIntegrationScheme = NULL;
  }

  ThermalSphericParticle::ThermalSphericParticle(IndexType NewId, NodesArrayType const& ThisNodes):SphericParticle(NewId, ThisNodes) {
    mpThermalIntegrationScheme = NULL;
  }

  ThermalSphericParticle::ThermalSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):SphericParticle(NewId, pGeometry, pProperties) {
    mpThermalIntegrationScheme = NULL;
  }

  Element::Pointer ThermalSphericParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
    return SphericParticle::Pointer(new ThermalSphericParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
  }

  ThermalSphericParticle::~ThermalSphericParticle() {
    if (mpThermalIntegrationScheme != NULL) {
      delete mpThermalIntegrationScheme;
      mpThermalIntegrationScheme = NULL;
    }
  }

  //=====================================================================================================================================================================================
  // Initialization methods

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::Initialize(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Initialize base class
    SphericParticle::Initialize(r_process_info);

    // Set thermal flags
    this->Set(DEMThermalFlags::HAS_MOTION,                       r_process_info[MOTION_OPTION]);
    this->Set(DEMThermalFlags::HAS_DIRECT_CONDUCTION,            r_process_info[DIRECT_CONDUCTION_OPTION]);
    this->Set(DEMThermalFlags::HAS_INDIRECT_CONDUCTION,          r_process_info[INDIRECT_CONDUCTION_OPTION]);
    this->Set(DEMThermalFlags::HAS_CONVECTION,                   r_process_info[CONVECTION_OPTION]);
    this->Set(DEMThermalFlags::HAS_RADIATION,                    r_process_info[RADIATION_OPTION]);
    this->Set(DEMThermalFlags::HAS_FRICTION_HEAT,                r_process_info[FRICTION_HEAT_OPTION]);
    this->Set(DEMThermalFlags::HAS_ADJUSTED_CONTACT,             r_process_info[ADJUSTED_CONTACT_OPTION]);
    this->Set(DEMThermalFlags::HAS_TEMPERATURE_DEPENDENT_RADIUS, r_process_info[TEMPERATURE_DEPENDENT_RADIUS_OPTION]);

    // Set thermal integration scheme
    ThermalDEMIntegrationScheme::Pointer& thermal_integration_scheme = GetProperties()[DEM_THERMAL_INTEGRATION_SCHEME_POINTER];
    SetThermalIntegrationScheme(thermal_integration_scheme);

    mStoreContactParam = this->Is(DEMThermalFlags::HAS_MOTION)        &&
                        (this->Is(DEMThermalFlags::HAS_FRICTION_HEAT) ||
                        (this->Is(DEMThermalFlags::HAS_DIRECT_CONDUCTION) && r_process_info[DIRECT_CONDUCTION_MODEL].compare("collisional") == 0));    

    // Clear maps
    mContactParamsParticle.clear();
    mContactParamsWall.clear();

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::InitializeSolutionStep(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Initialize base class
    if (this->Is(DEMThermalFlags::HAS_MOTION))
      SphericParticle::InitializeSolutionStep(r_process_info);

    // Check if it is time to evaluate thermal problem
    const int step = r_process_info[TIME_STEPS];
    const int freq = r_process_info[THERMAL_FREQUENCY];
    mIsTimeToSolve = (step > 0) && (freq != 0) && (step - 1) % freq == 0;

    // Number of steps passed between thermal evaluation steps
    mNumStepsEval = (r_process_info[TIME_STEPS] == 1) ? 1 : r_process_info[THERMAL_FREQUENCY];

    // Save pre-step temperature
    mPreviousTemperature = GetParticleTemperature();

    // Initialize number of contact particle neighbors
    // (currently used only for cleaning contact parameters map)
    mNumberOfContactParticleNeighbor = 0;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::InitializeHeatFluxComputation(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Initialize heat fluxes contributions
    mConductionDirectHeatFlux   = 0.0;
    mConductionIndirectHeatFlux = 0.0;
    mRadiationHeatFlux          = 0.0;
    mFrictionHeatFlux           = 0.0;
    mConvectionHeatFlux         = 0.0;
    mPrescribedHeatFlux         = 0.0;
    mTotalHeatFlux              = 0.0;

    // Initialize environment-related variables for radiation
    if (this->Is(DEMThermalFlags::HAS_RADIATION)) {
      mRadiativeNeighbors     = 0;
      mEnvironmentTemperature = 0.0;
      mEnvironmentTempAux     = 0.0;
    }

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Computation methods

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::CalculateRightHandSide(const ProcessInfo& r_process_info, double dt, const array_1d<double, 3>& gravity) {
    KRATOS_TRY

    // Force components
    if (this->Is(DEMThermalFlags::HAS_MOTION))
      SphericParticle::CalculateRightHandSide(r_process_info, dt, gravity);
    
    // Heat flux components
    if (mIsTimeToSolve)
      ComputeHeatFluxes(r_process_info);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputeHeatFluxes(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Initialize heat fluxes computation
    InitializeHeatFluxComputation(r_process_info);

    // Compute heat fluxes with neighbor particles
    for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
      if (mNeighbourElements[i] == NULL) continue;
      mNeighbor_p    = dynamic_cast<ThermalSphericParticle*>(mNeighbourElements[i]);
      mNeighborType  = PARTICLE_NEIGHBOR;
      mNeighborIndex = i;
      ComputeHeatFluxWithNeighbor(r_process_info);
    }

    // Compute heat fluxes with contact neighbor walls
    for (unsigned int i = 0; i < mNeighbourRigidFaces.size(); i++) {
      if (mNeighbourRigidFaces[i] == NULL) continue;
      mNeighbor_w    = dynamic_cast<DEMWall*>(mNeighbourRigidFaces[i]);
      mNeighborType  = WALL_NEIGHBOR_CONTACT;
      mNeighborIndex = i;
      ComputeHeatFluxWithNeighbor(r_process_info);
    }

    // Compute heat fluxes with noncontact neighbor walls
    // ATTENTION:
    // Only one element of a noncontact wall is elected to represent the entire wall,
    // which is the edge(2D)/face(3D) intersected by the normal vector between the wall and the particle.
    // This elected element represents an infinity wall.
    for (unsigned int i = 0; i < mNeighbourNonContactRigidFaces.size(); i++) {
      if (mNeighbourNonContactRigidFaces[i] == NULL) continue;
      mNeighbor_w = dynamic_cast<DEMWall*>(mNeighbourNonContactRigidFaces[i]);
      mNeighborType = WALL_NEIGHBOR_NONCONTACT;
      mNeighborIndex = i;
      ComputeHeatFluxWithNeighbor(r_process_info);
    }

    // Finalize radiation computation of continuous methods
    if (this->Is(DEMThermalFlags::HAS_RADIATION))
      ComputeContinuumRadiativeHeatFlux(r_process_info);

    // Compute convection with surrounding fluid
    if (this->Is(DEMThermalFlags::HAS_CONVECTION))
      ComputeConvectiveHeatFlux(r_process_info);

    // Prescribed heat flux
    ComputePrescribedHeatFlux(r_process_info);

    // Sum up heat fluxes contributions
    mTotalHeatFlux = mConductionDirectHeatFlux + mConductionIndirectHeatFlux + mRadiationHeatFlux + mFrictionHeatFlux + mConvectionHeatFlux + mPrescribedHeatFlux;
    SetParticleHeatFlux(mTotalHeatFlux);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::StoreBallToBallContactInfo(const ProcessInfo& r_process_info, SphericParticle::ParticleDataBuffer& data_buffer, SphericParticle* neighbor, double GlobalContactForce[3], bool sliding) {
    KRATOS_TRY

    if (!mStoreContactParam)
      return;

    // Increment number of contact particle neighbors
    mNumberOfContactParticleNeighbor++;

    // Local relavive velocity components (normal and tangential)
    std::vector<double> LocalRelativeVelocity{ data_buffer.mLocalRelVel[2], sqrt(data_buffer.mLocalRelVel[0] * data_buffer.mLocalRelVel[0] + data_buffer.mLocalRelVel[1] * data_buffer.mLocalRelVel[1]) };

    // Local contact force components (normal and sliding tangential)
    double LocalContactForce[3] = { 0.0 };
    GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, GlobalContactForce, LocalContactForce);
    std::vector<double> LocalForce{ LocalContactForce[2] };

    // Friction heat generation is not considered when particles are not sliding against each other, so tangent velocity is set to zero.
    // ATTENTION: Becareful when using the tangent velocity in other context that is not friction heat generation, as it can be zero.
    if (sliding)
      LocalForce.push_back(sqrt(LocalContactForce[0] * LocalContactForce[0] + LocalContactForce[1] * LocalContactForce[1]));
    else
      LocalForce.push_back(0.0);

    // Update contact parameters
    ContactParams params;
    params.updated_step   = r_process_info[TIME_STEPS];
    params.local_velocity = LocalRelativeVelocity;
    params.local_force    = LocalForce;

    // Keep impact parameters if contact is not new
    if (mContactParamsParticle.count(neighbor)) {
      params.impact_time     = mContactParamsParticle[neighbor].impact_time;
      params.impact_velocity = mContactParamsParticle[neighbor].impact_velocity;
    }
    // Set impact parameters for new contacts
    else {
      params.impact_time     = r_process_info[TIME];
      params.impact_velocity = LocalRelativeVelocity;
    }

    // Add/Update parameters in map
    mContactParamsParticle[neighbor] = params;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::StoreBallToRigidFaceContactInfo(const ProcessInfo& r_process_info, SphericParticle::ParticleDataBuffer& data_buffer, DEMWall* neighbor, double GlobalContactForce[3], bool sliding) {
    KRATOS_TRY

    if (!mStoreContactParam)
      return;

    // Local relavive velocity components (normal and tangential)
    std::vector<double> LocalRelativeVelocity{ data_buffer.mLocalRelVel[2], sqrt(data_buffer.mLocalRelVel[0] * data_buffer.mLocalRelVel[0] + data_buffer.mLocalRelVel[1] * data_buffer.mLocalRelVel[1]) };

    // Local contact force components (normal and sliding tangential)
    double LocalContactForce[3] = { 0.0 };
    GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, GlobalContactForce, LocalContactForce);
    std::vector<double> LocalForce{ LocalContactForce[2] };

    // Friction heat generation is not considered when particles are not sliding against each other, so tangent velocity is set to zero.
    // ATTENTION: Becareful when using the tangent velocity in other context that is not friction heat generation, as it can be zero.
    if (sliding)
      LocalForce.push_back(sqrt(LocalContactForce[0] * LocalContactForce[0] + LocalContactForce[1] * LocalContactForce[1]));
    else
      LocalForce.push_back(0.0);

    // Update contact parameters
    ContactParams params;
    params.updated_step   = r_process_info[TIME_STEPS];
    params.local_velocity = LocalRelativeVelocity;
    params.local_force    = LocalForce;

    // Keep impact parameters if contact is not new
    if (mContactParamsWall.count(neighbor)) {
      params.impact_time     = mContactParamsWall[neighbor].impact_time;
      params.impact_velocity = mContactParamsWall[neighbor].impact_velocity;
    }
    // Set impact parameters for new contacts
    else {
      params.impact_time     = r_process_info[TIME];
      params.impact_velocity = LocalRelativeVelocity;
    }

    // Add/Update parameters in map
    mContactParamsWall[neighbor] = params;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::Move(const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag) {
    // Time integration of motion
    if (this->Is(DEMThermalFlags::HAS_MOTION))
      SphericParticle::Move(delta_t, rotation_option, force_reduction_factor, StepFlag);

    // Time integration of temperature
    if (mIsTimeToSolve && !this->Is(DEMThermalFlags::HAS_FIXED_TEMPERATURE) && !this->Is(DEMThermalFlags::IS_ADIABATIC)) {
      GetThermalIntegrationScheme().UpdateTemperature(GetGeometry()[0], delta_t * mNumStepsEval, GetParticleHeatCapacity()); // TODO: Remove last argument (capacity becomes a node property accessed with GetFastProperties - same as MOMENT_OF_INERTIA)
    }
  }

  //=====================================================================================================================================================================================
  // Finalization methods

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::FinalizeSolutionStep(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    if (this->Is(DEMThermalFlags::HAS_MOTION))
      SphericParticle::FinalizeSolutionStep(r_process_info);

    // Remove non-contacting neighbors from maps of contact parameters
    if (mStoreContactParam)
      CleanContactParameters(r_process_info);

    // Update temperature/heat flux
    if (mIsTimeToSolve)
      UpdateTemperatureDependentRadius(r_process_info);

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Update methods

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::UpdateTemperatureDependentRadius(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    if (!this->Is(DEMThermalFlags::HAS_TEMPERATURE_DEPENDENT_RADIUS))
      return;

    // Update radius
    const double new_radius = GetParticleRadius() * (1.0 + GetParticleExpansionCoefficient() * (GetParticleTemperature() - mPreviousTemperature));
    SetParticleRadius(new_radius);

    // Update inertia
    SetParticleMomentInertia(CalculateMomentOfInertia());

    // TODO: update density

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion(const ProcessInfo& r_process_info,
                                                                                                double& thermalDeltDisp,
                                                                                                double& thermalRelVel,
                                                                                                SphericParticle* element2) {}

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(const ProcessInfo& r_process_info,
                                                                                              double DeltDisp[3], //IN GLOBAL AXES
                                                                                              double RelVel[3],   //IN GLOBAL AXES
                                                                                              double OldLocalCoordSystem[3][3],
                                                                                              double LocalCoordSystem[3][3],
                                                                                              SphericParticle* neighbor_iterator) {}

  //=====================================================================================================================================================================================
  // Heat fluxes computation

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputeHeatFluxWithNeighbor(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Check if neighbor is adiabatic
    if (CheckAdiabaticNeighbor())
      return;

    // Compute simulated or adjusted interaction properties
    ComputeInteractionProps(r_process_info);

    // Heat transfer mechanisms
    if (this->Is(DEMThermalFlags::HAS_DIRECT_CONDUCTION))   ComputeDirectConductionHeatFlux(r_process_info);
    if (this->Is(DEMThermalFlags::HAS_INDIRECT_CONDUCTION)) ComputeIndirectConductionHeatFlux(r_process_info);
    if (this->Is(DEMThermalFlags::HAS_RADIATION))           ComputeRadiativeHeatFlux(r_process_info);
    if (this->Is(DEMThermalFlags::HAS_FRICTION_HEAT))       ComputeFrictionHeatFlux(r_process_info);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputeInteractionProps(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Set simulated properties
    mNeighborDistance   = ComputeDistanceToNeighbor();
    mNeighborSeparation = ComputeSeparationToNeighbor();
    mNeighborInContact  = CheckSurfaceDistance(0.0);
    mContactRadius      = ComputeContactRadius();

    // Set adjusted contact properties
    if (!mNeighborInContact || !this->Is(DEMThermalFlags::HAS_ADJUSTED_CONTACT)) {
      mNeighborDistanceAdjusted   = mNeighborDistance;
      mNeighborSeparationAdjusted = mNeighborSeparation;
      mContactRadiusAdjusted      = mContactRadius;
    }
    else {
      // Compute adjusted contact radius according to selected model
      std::string model = r_process_info[ADJUSTED_CONTACT_MODEL];
      if      (model.compare("zhou")   == 0) mContactRadiusAdjusted = AdjustedContactRadiusZhou(r_process_info);
      else if (model.compare("lu")     == 0) mContactRadiusAdjusted = AdjustedContactRadiusLu(r_process_info);
      else if (model.compare("morris") == 0) mContactRadiusAdjusted = AdjustedContactRadiusMorris(r_process_info);

      // Compute adjusted distance/separation from adjusted contact radius
      mNeighborDistanceAdjusted   = ComputeDistanceToNeighborAdjusted();
      mNeighborSeparationAdjusted = ComputeSeparationToNeighborAdjusted();
    }

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputeDirectConductionHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Check for contact
    if (!mNeighborInContact)
      return;

    // Compute heat flux according to selected model
    std::string model = r_process_info[DIRECT_CONDUCTION_MODEL];

    if      (model.compare("batchelor_obrien") == 0) mConductionDirectHeatFlux += DirectConductionBatchelorOBrien(r_process_info);
    else if (model.compare("thermal_pipe")     == 0) mConductionDirectHeatFlux += DirectConductionThermalPipe(r_process_info);
    else if (model.compare("collisional")      == 0) mConductionDirectHeatFlux += DirectConductionCollisional(r_process_info);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputeIndirectConductionHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Compute heat flux according to selected model
    std::string model = r_process_info[INDIRECT_CONDUCTION_MODEL];

    if      (model.compare("surrounding_layer") == 0) mConductionIndirectHeatFlux += IndirectConductionSurroundingLayer(r_process_info);
    else if (model.compare("voronoi_a")         == 0) mConductionIndirectHeatFlux += IndirectConductionVoronoiA(r_process_info);
    else if (model.compare("voronoi_b")         == 0) mConductionIndirectHeatFlux += IndirectConductionVoronoiB(r_process_info);
    else if (model.compare("vargas_mccarthy")   == 0) mConductionIndirectHeatFlux += IndirectConductionVargasMcCarthy(r_process_info);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputeRadiativeHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // TODO: radiation with walls not yet implemented
    if (mNeighborType & WALL_NEIGHBOR)
      return;

    // Check if particles are close enough
    if (!CheckSurfaceDistance(r_process_info[MAX_RADIATION_DISTANCE]))
      return;

    // Update number of radiative neighbors
    mRadiativeNeighbors++;

    // Accumulate environment temperature (in continuum methods) or compute heat flux (in discrete methods) according to selected model
    std::string model = r_process_info[RADIATION_MODEL];

    if (model.compare("continuum_zhou") == 0) {
      mEnvironmentTemperature += GetNeighborTemperature();
    }
    else if (model.compare("continuum_krause") == 0) {
      const double neighbor_emissivity  = GetNeighborEmissivity();
      const double neighbor_temperature = GetNeighborTemperature();
      const double neighbor_surface     = mNeighbor_p->GetParticleSurfaceArea();
      mEnvironmentTemperature += 0.5 * STEFAN_BOLTZMANN * neighbor_emissivity * neighbor_surface * pow(neighbor_temperature, 4.0);
      mEnvironmentTempAux     += 0.5 * STEFAN_BOLTZMANN * neighbor_emissivity * neighbor_surface;
    }

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputeContinuumRadiativeHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Check if radiation neighbors exist
    if (mRadiativeNeighbors == 0)
      return;

    // compute heat flux of continuous methods according to selected model
    std::string model = r_process_info[RADIATION_MODEL];

    if      (model.compare("continuum_zhou")   == 0) mRadiationHeatFlux += RadiationContinuumZhou(r_process_info);
    else if (model.compare("continuum_krause") == 0) mRadiationHeatFlux += RadiationContinuumKrause(r_process_info);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputeFrictionHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Check for contact
    if (!mNeighborInContact || !this->Is(DEMThermalFlags::HAS_MOTION))
      return;

    // Compute heat generation by friction according to selected model
    mFrictionHeatFlux += FrictionGenerationSlidingVelocity(r_process_info);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputeConvectiveHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY
    
    const double surface_area       = GetParticleSurfaceArea();
    const double char_length        = GetParticleCharacteristicLength();
    const double fluid_conductivity = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double temp_grad          = r_process_info[FLUID_TEMPERATURE] - GetParticleTemperature();

    // Compute Nusselt number according to selected model
    double Nu = 0.0;
    std::string model = r_process_info[CONVECTION_MODEL];

    if      (model.compare("sphere_hanz_marshall") == 0) Nu = NusseltHanzMarshall(r_process_info);
    else if (model.compare("sphere_whitaker")      == 0) Nu = NusseltWhitaker(r_process_info);
    else if (model.compare("sphere_gunn")          == 0) Nu = NusseltGunn(r_process_info);
    else if (model.compare("sphere_li_mason")      == 0) Nu = NusseltLiMason(r_process_info);

    // Compute heat flux
    mConvectionHeatFlux += (Nu * fluid_conductivity / char_length) * surface_area * temp_grad;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputePrescribedHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY
    
    // Heat flux over surface area
    if (mPrescribedHeatFluxSurface != 0.0)
      mPrescribedHeatFlux += mPrescribedHeatFluxSurface * GetParticleSurfaceArea();

    // Volume heat source
    if (mPrescribedHeatFluxVolume != 0.0)
      mPrescribedHeatFlux += mPrescribedHeatFluxVolume * GetParticleVolume();

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Heat transfer models

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::DirectConductionBatchelorOBrien(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    const double keff      = ComputeEffectiveConductivity();
    const double temp_grad = GetNeighborTemperature() - GetParticleTemperature();

    return 4.0 * keff * mContactRadiusAdjusted * temp_grad;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::DirectConductionThermalPipe(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    const double kavg      = ComputeAverageConductivity();
    const double temp_grad = GetNeighborTemperature() - GetParticleTemperature();

    return kavg * (Globals::Pi * mContactRadiusAdjusted * mContactRadiusAdjusted) * temp_grad / mNeighborDistanceAdjusted;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::DirectConductionCollisional(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Get collision time and impact normal velocity
    typename ContactParams contact_params = GetContactParameters();
    const double col_time                 = r_process_info[TIME] - contact_params.impact_time;
    const double impact_normal_velocity   = fabs(contact_params.impact_velocity[0]);

    // Compute max collision time
    double col_time_max = 0.0;
    if (impact_normal_velocity != 0.0)
      col_time_max = ComputeMaxCollisionTime();
    
    // Check if collision time is smaller than max value, otherwise use static model (batchelor_obrien)
    if (col_time < col_time_max) {
      const double temp_grad = GetNeighborTemperature() - GetParticleTemperature();
      const double Rc_max    = ComputeMaxContactRadius(); // TODO: This should be multiplied by the correction coefficient (and not computed with real Young modulus)
      const double Fo        = ComputeFourierNumber();

      const double a1 = GetParticleDensity() * GetParticleHeatCapacity();
      const double a2 = GetNeighborDensity() * GetNeighborHeatCapacity();
      const double b1 = a1 * GetParticleConductivity();
      const double b2 = a2 * GetNeighborConductivity();
      const double c  = a1 / a2;

      const double C1 = -2.300 * c * c +  8.909 * c - 4.235;
      const double C2 =  8.169 * c * c - 33.770 * c + 24.885;
      const double C3 = -5.758 * c * c + 24.464 * c - 20.511;

      const double C_coeff = 0.435 * (sqrt(C2 * C2 - 4.0 * C1 * (C3 - Fo)) - C2) / C1;

      return C_coeff * Globals::Pi * Rc_max * Rc_max * pow(col_time_max,-0.5) * temp_grad / (pow(b1,-0.5) + pow(b2,-0.5));
    }
    else {
      return DirectConductionBatchelorOBrien(r_process_info);
    }

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::IndirectConductionSurroundingLayer(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Check if particles are close enough
    const double layer = r_process_info[FLUID_LAYER_THICKNESS];
    if (!CheckSurfaceDistance(layer))
      return 0.0;

    // Compute heat transfer coefficient
    const double fluid_conductivity = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double min_dist           = r_process_info[MIN_CONDUCTION_DISTANCE];
    double h = 0.0;

    if (mNeighborType & PARTICLE_NEIGHBOR) {
      const double particle_radius = GetParticleRadius();
      const double neighbor_radius = GetNeighborRadius();
      const double r_min           = std::min(particle_radius, neighbor_radius);
      const double r_max           = std::max(particle_radius, neighbor_radius);

      // Compute upper limit of integral
      const double param = pow((r_max + (layer * r_max)), 2.0);
      double upp_lim;

      if (mNeighborDistanceAdjusted <= sqrt(param - r_min * r_min))
        upp_lim = r_min;
      else
        upp_lim = sqrt(param - pow(((param - r_min * r_min + mNeighborDistanceAdjusted * mNeighborDistanceAdjusted) / (2.0 * mNeighborDistanceAdjusted)), 2.0));

      // Build struct of integration parameters
      IntegrandParams params;
      params.p1 = min_dist;
      params.p2 = particle_radius;
      params.p3 = neighbor_radius;

      // Heat transfer coefficient from integral expression solved numerically
      h = fluid_conductivity * AdaptiveSimpsonIntegration(r_process_info, mContactRadiusAdjusted, upp_lim, params, &ThermalSphericParticle::EvalIntegrandSurrLayer);
    }
    else if (mNeighborType & WALL_NEIGHBOR) {
      const double particle_radius = GetParticleRadius();
      double a, b, c, r_in, r_out;

      a = (mNeighborDistanceAdjusted - particle_radius) / particle_radius;

      if (mNeighborDistance > particle_radius + min_dist)
        r_in = 0.0;
      else
        r_in = sqrt(1.0 - pow(min_dist / particle_radius - a - 1.0, 2.0));

      if (a > sqrt(pow((particle_radius + (layer * particle_radius)) / particle_radius, 2.0) - 1.0) - 1.0)
        r_out = sqrt(pow((particle_radius + (layer * particle_radius)) / particle_radius, 2.0) - pow(a + 1.0, 2.0));
      else
        r_out = 1.0;

      b = sqrt(1.0 - r_out * r_out);
      c = sqrt(1.0 - r_in  * r_in);

      // Heat transfer coefficient from analytical solution of the integral expression
      h = 2.0 * Globals::Pi * fluid_conductivity * particle_radius * ((a + 1.0) * log(fabs(b - a - 1.0) / fabs(a - c + 1.0)) + b - c);
    }

    // Compute heat flux
    return h * (GetNeighborTemperature() - GetParticleTemperature());

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::IndirectConductionVoronoiA(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Check if particles are close enough
    if (!CheckSurfaceDistance(r_process_info[MAX_CONDUCTION_DISTANCE]))
      return 0.0;

    // Get radii
    const double particle_radius = GetParticleRadius();
    const double neighbor_radius = GetNeighborRadius();

    // Get radius of voronoi cell face
    double rij = GetVoronoiCellFaceRadius(r_process_info);
    if (rij <= mContactRadiusAdjusted)
      return 0.0;

    // Compute heat transfer coefficient
    double h = 0.0;

    if (mNeighborType & WALL_NEIGHBOR) {
      // Compute upper limit of integral
      double upp_lim = particle_radius * rij / sqrt(rij * rij + mNeighborDistanceAdjusted * mNeighborDistanceAdjusted);

      // Build struct of integration parameters
      IntegrandParams params;
      params.p1 = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
      params.p2 = GetParticleConductivity();
      params.p3 = particle_radius;
      params.p4 = rij;

      // Heat transfer coefficient from integral expression solved numerically
      h = AdaptiveSimpsonIntegration(r_process_info, mContactRadiusAdjusted, upp_lim, params, &ThermalSphericParticle::EvalIntegrandVoronoiWall);
    }
    else if (particle_radius == neighbor_radius) {
      // Compute upper limit of integral
      double upp_lim = particle_radius * rij / sqrt(rij * rij + mNeighborDistanceAdjusted * mNeighborDistanceAdjusted / 4.0);

      // Build struct of integration parameters
      IntegrandParams params;
      params.p1 = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
      params.p2 = ComputeEffectiveConductivity();
      params.p3 = particle_radius;
      params.p4 = rij;

      // Heat transfer coefficient from integral expression solved numerically
      h = AdaptiveSimpsonIntegration(r_process_info, mContactRadiusAdjusted, upp_lim, params, &ThermalSphericParticle::EvalIntegrandVoronoiMono);
    }
    else {
      double D1, D2, rij_, upp_lim;

      if (mNeighborInContact)
        D1 = sqrt(particle_radius * particle_radius - mContactRadiusAdjusted * mContactRadiusAdjusted);
      else
        D1 = (particle_radius * particle_radius - neighbor_radius * neighbor_radius + mNeighborDistanceAdjusted * mNeighborDistanceAdjusted) / (2 * mNeighborDistanceAdjusted);

      D2 = mNeighborDistanceAdjusted - D1;

      if (particle_radius <= neighbor_radius)
        upp_lim = particle_radius * rij / sqrt(rij * rij + D1 * D1);
      else
        upp_lim = neighbor_radius * rij / sqrt(rij * rij + D2 * D2);

      rij_ = D2 * upp_lim / sqrt(neighbor_radius * neighbor_radius - upp_lim * upp_lim);

      // Build struct of integration parameters
      IntegrandParams params;
      params.p1 = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
      params.p2 = GetParticleConductivity();
      params.p3 = GetNeighborConductivity();
      params.p4 = particle_radius;
      params.p5 = neighbor_radius;
      params.p6 = rij;
      params.p7 = rij_;
      params.p8 = D1;
      params.p9 = D2;

      // Heat transfer coefficient from integral expression solved numerically
      h = AdaptiveSimpsonIntegration(r_process_info, mContactRadiusAdjusted, upp_lim, params, &ThermalSphericParticle::EvalIntegrandVoronoiMulti);
    }

    // Compute heat flux
    return h * (GetNeighborTemperature() - GetParticleTemperature());

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::IndirectConductionVoronoiB(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Check if particles are close enough
    if (!CheckSurfaceDistance(r_process_info[MAX_CONDUCTION_DISTANCE]))
      return 0.0;

    // Get parameters
    const double particle_radius       = GetParticleRadius();
    const double neighbor_radius       = GetNeighborRadius();
    const double particle_conductivity = GetParticleConductivity();
    const double neighbor_conductivity = GetNeighborConductivity();
    const double fluid_conductivity    = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double core                  = r_process_info[ISOTHERMAL_CORE_RADIUS];

    // Get radius of voronoi cell face
    const double rij = GetVoronoiCellFaceRadius(r_process_info);
    if (rij <= mContactRadiusAdjusted)
      return 0.0;

    // Compute heat transfer coefficient
    double h = 0.0;

    if (mNeighborType & WALL_NEIGHBOR) {
      const double kp = GetParticleConductivity();
      const double rc = core * particle_radius;
      const double d  = mNeighborDistanceAdjusted;
      const double a  = (1.0 / rc - 1.0 / particle_radius) / (2.0 * kp) + 1.0 / (2 * fluid_conductivity * particle_radius);
      const double b  = 1.0 / (2 * fluid_conductivity * d);
      const double c0 = d / sqrt(rij * rij + d * d);
      const double c1 = d / sqrt(mContactRadiusAdjusted * mContactRadiusAdjusted + d * d);
      const double f  = (a - b * c0) / (a - b * c1);
      double ln = 0.0;
      if (f > 0.0)
        ln = log(f);

      // Heat transfer coefficient
      h = Globals::Pi * ln / b;
    }
    else if (particle_radius == neighbor_radius) {
      const double keff = ComputeEffectiveConductivity();
      const double rc   = core * particle_radius;
      const double D    = mNeighborDistanceAdjusted / 2.0;
      const double a    = (1.0 / rc - 1.0 / particle_radius) / (2.0 * keff) + 1.0 / (fluid_conductivity * particle_radius);
      const double b    = 1.0 / (fluid_conductivity * D);
      const double c0   = D / sqrt(rij * rij + D * D);
      const double c1   = D / sqrt(mContactRadiusAdjusted * mContactRadiusAdjusted + D * D);
      const double f    = (a - b * c0) / (a - b * c1);
      double ln = 0.0;
      if (f > 0.0)
        ln = log(f);

      // Heat transfer coefficient
      h = Globals::Pi * ln / b;
    }
    else {
      const double An = Globals::Pi * rij * rij; // area of neighboring voronoi cells

      const double gamma1 = particle_radius / mNeighborDistanceAdjusted;
      const double gamma2 = neighbor_radius / mNeighborDistanceAdjusted;
      const double dgamma = gamma2 - gamma1;

      const double A = (particle_conductivity + fluid_conductivity * (1.0 / core - 1.0)) / (particle_conductivity * gamma1);
      const double B = (neighbor_conductivity + fluid_conductivity * (1.0 / core - 1.0)) / (neighbor_conductivity * gamma2);

      const double lambda = (1.0 + dgamma * A) * (1.0 - dgamma * B);

      const double delmax = 0.5 * (sqrt((4.0 * An) / (Globals::Pi * mNeighborDistanceAdjusted * mNeighborDistanceAdjusted * (1.0 - dgamma * dgamma)) + 1.0) - dgamma);
      const double delmin = 0.5 * (sqrt((4.0 * mContactRadiusAdjusted * mContactRadiusAdjusted) / (mNeighborDistanceAdjusted * mNeighborDistanceAdjusted * (1.0 - dgamma * dgamma)) + 1.0) - dgamma);

      const double Xmax = ((A + B) * delmax + dgamma * B - 1.0) / sqrt(fabs(lambda));
      const double Xmin = ((A + B) * delmin + dgamma * B - 1.0) / sqrt(fabs(lambda));

      const double Y1 = (Xmax - Xmin) / (1.0 - Xmax * Xmin);
      const double Y2 = (Xmax - Xmin) / (1.0 + Xmax * Xmin);

      // Heat transfer coefficient
      if (lambda > 0.0)
        h = Globals::Pi * fluid_conductivity * mNeighborDistanceAdjusted * (1.0 - dgamma * dgamma) * log(fabs((1.0 - Y1) / (1.0 + Y1))) / (2.0 * sqrt(fabs(lambda)));
      else if (lambda < 0.0)
        h = Globals::Pi * fluid_conductivity * mNeighborDistanceAdjusted * (1.0 - dgamma * dgamma) * atan(Y2) / (2.0 * sqrt(fabs(lambda)));
      else
        h = Globals::Pi * fluid_conductivity * mNeighborDistanceAdjusted * (1.0 - dgamma * dgamma) * (1.0 / delmin - 1.0 / delmax) / (A + B);
    }

    // Compute heat flux
    return h * (GetNeighborTemperature() - GetParticleTemperature());

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::IndirectConductionVargasMcCarthy(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Check for contact
    if (!mNeighborInContact)
      return 0.0;

    // Assumption 1: Formulation for a liquid (not gas) as the interstitial fluid is being used
    const double fluid_conductivity = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double temp_grad          = GetNeighborTemperature() - GetParticleTemperature();

    // Assumption 2 : Model developed for mono-sized particles, but the average radius is being used (if neighbor is a wall, it is assumed as a particle with the same radius)
    const double particle_radius = GetParticleRadius();
    const double neighbor_radius = (mNeighborType & PARTICLE_NEIGHBOR) ? GetNeighborRadius() : particle_radius;
    const double avg_radius      = (particle_radius + neighbor_radius) / 2.0;

    // Compute heat flux
    return 4.0 * Globals::Pi * fluid_conductivity * (1.0 - 0.5 * pow(mContactRadiusAdjusted / avg_radius, 2.0) * (avg_radius - mContactRadiusAdjusted)) * temp_grad / (1.0 - Globals::Pi / 4.0);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::NusseltHanzMarshall(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    const double Pr = ComputePrandtlNumber(r_process_info);
    const double Re = ComputeReynoldNumber(r_process_info);

    return 2.0 + 0.6 * pow(Re,0.5) * pow(Pr,1.0/3.0);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::NusseltWhitaker(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    const double Pr = ComputePrandtlNumber(r_process_info);
    const double Re = ComputeReynoldNumber(r_process_info);

    // Assumption: temperature-dependent viscosity at particle surface is negleted
    return 2.0 + (0.4 * pow(Re,0.5) + 0.06 * pow(Re,2.0/3.0)) * pow(Pr,0.4);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::NusseltGunn(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    const double Pr  = ComputePrandtlNumber(r_process_info);
    const double Re  = ComputeReynoldNumber(r_process_info);
    const double por = r_process_info[AVERAGE_POROSITY];
    
    return (7.0 - 10.0 * por + 5.0 * por * por) * (1.0 + 0.7 * pow(Re,0.2) * pow(Pr,1.0/3.0)) + (1.33 - 2.4 * por + 1.2 * por * por) * pow(Re,0.7) * pow(Pr,1.0/3.0);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::NusseltLiMason(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    const double Pr  = ComputePrandtlNumber(r_process_info);
    const double Re  = ComputeReynoldNumber(r_process_info);
    const double por = r_process_info[AVERAGE_POROSITY];
    const double m   = 4.75; // Assumption: exponent "m = 4.75" recommended for dense systems (3.50 is recommended for dilute systems)

    if      (Re < 200.0)  return 2.0 + 0.6 * pow(por,m) * pow(Re,0.5) * pow(Pr,1.0/3.0);
    else if (Re < 1500.0) return 2.0 + pow(por,m) * (0.5 * pow(Re,0.5) + 0.02 * pow(Re,0.8)) * pow(Pr,1.0/3.0);
    else                  return 2.0 + 0.000045 * pow(por,m) * pow(Re,1.8);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::RadiationContinuumZhou(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Get parameters
    const double particle_emissivity  = GetParticleEmissivity();
    const double particle_surface     = GetParticleSurfaceArea();
    const double particle_temperature = GetParticleTemperature();
    const double porosity             = r_process_info[AVERAGE_POROSITY];
    const double f_temperature        = r_process_info[FLUID_TEMPERATURE];

    // Compute final value of environment temperature
    const double env_temperature = porosity * f_temperature + (1.0 - porosity) * mEnvironmentTemperature / mRadiativeNeighbors;

    // Compute heat flux
    return STEFAN_BOLTZMANN * particle_emissivity * particle_surface * (pow(env_temperature,4.0) - pow(particle_temperature,4.0));

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::RadiationContinuumKrause(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Get parameters
    const double particle_emissivity  = GetParticleEmissivity();
    const double particle_surface     = GetParticleSurfaceArea();
    const double particle_temperature = GetParticleTemperature();

    // Compute final value of environment temperature
    const double env_temperature = pow(mEnvironmentTemperature / mEnvironmentTempAux, 0.25);

    // Compute heat flux
    return STEFAN_BOLTZMANN * particle_emissivity * particle_surface * (pow(env_temperature,4.0) - pow(particle_temperature,4.0));

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::FrictionGenerationSlidingVelocity(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Model parameters
    typename ContactParams contact_params = GetContactParameters();
    const double velocity_tangent         = contact_params.local_velocity[1];
    const double force_normal             = contact_params.local_force[0];
    if (velocity_tangent == 0 || force_normal == 0) return 0.0;
    const double friction_conversion      = r_process_info[FRICTION_HEAT_CONVERSION];
    const double friction_coeff           = GetContactDynamicFrictionCoefficient();

    // Partition coefficient
    const double k1 = GetParticleConductivity();
    const double k2 = GetNeighborConductivity();
    const double partition = k1 / (k1 + k2);
    
    // Compute frictional heat transfer
    return partition * friction_conversion * friction_coeff * fabs(velocity_tangent * force_normal);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::AdjustedContactRadiusZhou(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Simulation and real values of effective Young modulus
    const double eff_young      = ComputeEffectiveYoung();
    const double eff_young_real = ComputeEffectiveYoungReal();

    // Adjusted value of contact radius
    return mContactRadius * pow(eff_young / eff_young_real, 0.2);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::AdjustedContactRadiusLu(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Effective radius
    const double eff_radius = ComputeEffectiveRadius();

    // Simulation and real values of effective Young modulus
    const double eff_young      = ComputeEffectiveYoung();
    const double eff_young_real = ComputeEffectiveYoungReal();

    // Simulation and real values of stiffness
    const double stiff      = 4.0 / 3.0 * sqrt(eff_radius) * eff_young;
    const double stiff_real = 4.0 / 3.0 * sqrt(eff_radius) * eff_young_real;

    // Adjusted value of contact radius
    return pow(mContactRadius * stiff / stiff_real, 2.0/3.0);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::AdjustedContactRadiusMorris(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Parameters
    const double eff_young      = ComputeEffectiveYoung();
    const double eff_young_real = ComputeEffectiveYoungReal();
    const double eff_radius     = ComputeEffectiveRadius();
    const double identation     = std::max(mNeighborSeparation, 0.0);

    // Contact force with simulation parameters (using Hertz theory)
    const double hertz_force = 4.0 * eff_young * sqrt(eff_radius) * pow(identation, 3.0 / 2.0) / 3.0;

    // Area correction
    const double correction_area = pow(hertz_force * eff_radius / eff_young_real, 1.0 / 3.0);

    // Time correction
    double correction_time = 1.0;

    if (this->Is(DEMThermalFlags::HAS_MOTION)) {
      // TODO: Compute time correction from collision time
      correction_time = 1.0;
    }

    // Adjusted value of contact radius
    return correction_area * correction_time;

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Auxiliary computations

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputeAddedSearchDistance(const ProcessInfo& r_process_info, double& added_search_distance) {
    KRATOS_TRY

    if (this->Is(DEMThermalFlags::HAS_INDIRECT_CONDUCTION)) {
      std::string model = r_process_info[INDIRECT_CONDUCTION_MODEL];
      if (model.compare("surrounding_layer") == 0) {
        const double model_search_distance  = GetParticleRadius() * r_process_info[FLUID_LAYER_THICKNESS];
        const double current_added_distance = added_search_distance;
        added_search_distance = std::max(current_added_distance, model_search_distance);
      }
      else if (model.compare("voronoi_a") == 0 ||
               model.compare("voronoi_b") == 0) {
        const double model_search_distance  = GetParticleRadius() * r_process_info[MAX_CONDUCTION_DISTANCE];
        const double current_added_distance = added_search_distance;
        added_search_distance = std::max(current_added_distance, model_search_distance);
      }
    }

    if (this->Is(DEMThermalFlags::HAS_RADIATION)) {
      std::string model = r_process_info[RADIATION_MODEL];
      if (model.compare("continuum_zhou")   == 0 ||
          model.compare("continuum_krause") == 0) {
        const double model_search_distance  = GetParticleRadius() * (r_process_info[MAX_RADIATION_DISTANCE]);
        const double current_added_distance = added_search_distance;
        added_search_distance = std::max(current_added_distance, model_search_distance);
      }
    }

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputePrandtlNumber(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    const double fluid_viscosity    = r_process_info[FLUID_VISCOSITY];
    const double fluid_heatcapacity = r_process_info[FLUID_HEAT_CAPACITY];
    const double fluid_conductivity = r_process_info[FLUID_THERMAL_CONDUCTIVITY];

    return fluid_viscosity * fluid_heatcapacity / fluid_conductivity;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeReynoldNumber(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    const double char_length     = GetParticleCharacteristicLength();
    const double rel_velocity    = ComputeFluidRelativeVelocity(r_process_info);
    const double fluid_density   = r_process_info[FLUID_DENSITY];
    const double fluid_viscosity = r_process_info[FLUID_VISCOSITY];

    return fluid_density * char_length * rel_velocity / fluid_viscosity;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeFluidRelativeVelocity(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    array_1d<double, 3> rel_velocity;
    noalias(rel_velocity) = GetParticleVelocity() - r_process_info[FLUID_VELOCITY];
    return DEM_MODULUS_3(rel_velocity);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetVoronoiCellFaceRadius(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    const double particle_radius = GetParticleRadius();
    const double neighbor_radius = GetNeighborRadius();

    // Based on voronoi diagram
    if (r_process_info[VORONOI_METHOD].compare("tesselation") == 0) {
      if (mNeighborType & WALL_NEIGHBOR) {
        // Assumption: radius of voronoi cell face is proportional to the particle radius
        return 1.5 * particle_radius;
      }
      else if (mNeighborVoronoiRadius.count(mNeighbor_p->mDelaunayPointListIndex)) {
        return mNeighborVoronoiRadius[mNeighbor_p->mDelaunayPointListIndex];
      }
      else {
        return 0.0;
      }
    }

    // Based on porosity
    else if (r_process_info[VORONOI_METHOD].compare("posority") == 0) {
      if (mNeighborType & WALL_NEIGHBOR) {
        // Assumption: using particle radius only
        return 0.56 * particle_radius * pow((1.0 - r_process_info[AVERAGE_POROSITY]), -1.0 / 3.0);
      }
      else {
        // Assumption: using average radius
        return 0.56 * (particle_radius + neighbor_radius) / 2.0 * pow((1.0 - r_process_info[AVERAGE_POROSITY]), -1.0 / 3.0);
      }
    }

    else {
      return 0.0;
    }

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Numerical integration

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::AdaptiveSimpsonIntegration(const ProcessInfo& r_process_info, double a, double b, IntegrandParams params, double (ThermalSphericParticle::*evalIntegrand)(IntegrandParams)) {
    KRATOS_TRY

    // Initialization
    params.x = a;
    double fa = (this->*evalIntegrand)(params);

    params.x = b;
    double fb = (this->*evalIntegrand)(params);

    params.x = (a + b) / 2.0;
    double fc = (this->*evalIntegrand)(params);

    // Get tolerance
    double tol = r_process_info[INTEGRAL_TOLERANCE];
    constexpr double eps = std::numeric_limits<double>::epsilon();
    if (tol < 10.0 * eps) tol = 10.0 * eps;

    // Solve integral recursively with adaptive Simpson quadrature
    return RecursiveSimpsonIntegration(a, b, fa, fb, fc, tol, params, evalIntegrand);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::RecursiveSimpsonIntegration(double a, double b, double fa, double fb, double fc, double tol, IntegrandParams params, double (ThermalSphericParticle::*evalIntegrand)(IntegrandParams)) {
    KRATOS_TRY

    // TODO: in order to catch possible erros that can occur in singularities,
    //       add a min value for subdivision size (to contain machine representable point) and a max number of function evaluation (+- 10000).

    double c = (a + b) / 2.0;

    params.x = (a + c) / 2.0;
    double fd = (this->*evalIntegrand)(params);

    params.x = (c + b) / 2.0;
    double fe = (this->*evalIntegrand)(params);

    double I1 = (b - a) / 6.0  * (fa + 4.0 * fc + fb);
    double I2 = (b - a) / 12.0 * (fa + 4.0 * fd + 2.0 * fc + 4.0 * fe + fb);

    if (fabs(I2 - I1) <= tol) {
      return I2 + (I2 - I1) / 15.0;
    }
    else { // sub-divide interval recursively
      double Ia = RecursiveSimpsonIntegration(a, c, fa, fc, fd, tol, params, evalIntegrand);
      double Ib = RecursiveSimpsonIntegration(c, b, fc, fb, fe, tol, params, evalIntegrand);
      return Ia + Ib;
    }

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::EvalIntegrandSurrLayer(IntegrandParams params) {
    KRATOS_TRY

    const double r    = params.x;
    const double dmin = params.p1;
    const double r1   = params.p2;
    const double r2   = params.p3;

    return 2.0 * Globals::Pi * r / std::max(dmin, mNeighborDistanceAdjusted - sqrt(r1 * r1 - r * r) - sqrt(r2 * r2 - r * r));

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::EvalIntegrandVoronoiWall(IntegrandParams params) {
    KRATOS_TRY

    const double r    = params.x;
    const double kf   = params.p1;
    const double kp   = params.p2;
    const double rp   = params.p3;
    const double rij  = params.p4;

    return 2.0 * Globals::Pi * r / ((sqrt(rp * rp - r * r) - r * mNeighborDistanceAdjusted / rij) / kp + (mNeighborDistanceAdjusted - sqrt(rp * rp - r * r)) / kf);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::EvalIntegrandVoronoiMono(IntegrandParams params) {
    KRATOS_TRY

    const double r    = params.x;
    const double kf   = params.p1;
    const double keff = params.p2;
    const double rp   = params.p3;
    const double rij  = params.p4;

    return 2.0 * Globals::Pi * r / ((sqrt(rp * rp - r * r) - r * mNeighborDistanceAdjusted / (2.0 * rij)) / keff + 2.0 * (mNeighborDistanceAdjusted / 2.0 - sqrt(rp * rp - r * r)) / kf);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::EvalIntegrandVoronoiMulti(IntegrandParams params) {
    KRATOS_TRY

    const double r    = params.x;
    const double kf   = params.p1;
    const double k1   = params.p2;
    const double k2   = params.p3;
    const double r1   = params.p4;
    const double r2   = params.p5;
    const double rij  = params.p6;
    const double rij_ = params.p7;
    const double D1   = params.p8;
    const double D2   = params.p9;

    const double beta1 = sqrt(r1 * r1 - r * r);
    const double beta2 = sqrt(r2 * r2 - r * r);

    return 2.0 * Globals::Pi * r / ((beta1 - D1 * r / rij) / k1 + (beta2 - D2 * r / rij_) / k2 + (mNeighborDistanceAdjusted - beta1 - beta2) / kf);

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Neighbor interaction computations

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::CleanContactParameters(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // When size of contact parameters map is different from the number of contacting neighbors,
    // it means that there are additional neighbors in the map that are no longer in contact and must be removed.
    // PS1: mNumberOfContactParticleNeighbor must count the number of contacting particle neighbors.
    // PS2: An indication that a neighbor must be removed from the map is when the updated_step parameter is not the current.

    if (mContactParamsParticle.size() != mNumberOfContactParticleNeighbor)
      for (auto it = mContactParamsParticle.cbegin(); it != mContactParamsParticle.cend(); )
        it = (it->second.updated_step != r_process_info[TIME_STEPS]) ? mContactParamsParticle.erase(it) : std::next(it);

    if (mContactParamsWall.size() != mNeighbourRigidFaces.size())
      for (auto it = mContactParamsWall.cbegin(); it != mContactParamsWall.cend(); )
        it = (it->second.updated_step != r_process_info[TIME_STEPS]) ? mContactParamsWall.erase(it) : std::next(it);
    
    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  bool ThermalSphericParticle::CheckAdiabaticNeighbor() {
    KRATOS_TRY

    return ((mNeighborType & PARTICLE_NEIGHBOR && mNeighbor_p->Is(DEMThermalFlags::IS_ADIABATIC)) ||
            (mNeighborType & WALL_NEIGHBOR     && mNeighbor_w->Is(DEMThermalFlags::IS_ADIABATIC)));

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  bool ThermalSphericParticle::CheckSurfaceDistance(const double radius_factor) {
    KRATOS_TRY

    // Assumption: radius_factor applies to the larger radius
    const double particle_radius = GetParticleRadius();
    const double neighbor_radius = GetNeighborRadius(); // must be zero for walls
    return (mNeighborSeparation < radius_factor * std::max(particle_radius, neighbor_radius));

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeDistanceToNeighbor() {
    KRATOS_TRY

    if (mNeighborType & PARTICLE_NEIGHBOR) {
      array_1d<double, 3> direction;
      noalias(direction) = GetParticleCoordinates() - GetNeighborCoordinates();
      return DEM_MODULUS_3(direction);
    }
    else if (mNeighborType & WALL_NEIGHBOR_CONTACT) {
      // Computing the distance again, as it is done in SphericParticle::ComputeBallToRigidFaceContactForce
      double distance = 0.0;
      array_1d<double, 4>& weight = this->mContactConditionWeights[mNeighborIndex];

      // Dummy variables: not used now
      double dummy1[3][3];
      DEM_SET_COMPONENTS_TO_ZERO_3x3(dummy1);
      array_1d<double, 3> dummy2 = ZeroVector(3);
      array_1d<double, 3> dummy3 = ZeroVector(3);
      int dummy4 = 0;

      mNeighbor_w->ComputeConditionRelativeData(mNeighborIndex, this, dummy1, distance, weight, dummy2, dummy3, dummy4);

      return distance;
    }
    else if (mNeighborType & WALL_NEIGHBOR_NONCONTACT) {
      // Computing the distance again, as it is done in SphericParticle::ComputeBallToRigidFaceContactForce
      double distance = 0.0;

      // ATTENTION:
      // Weight vector defines the indexes of the wall nodes in RigidEdge2D::ComputeConditionRelativeData.
      // It is also used to compute the number of points in RigidFace3D::ComputeConditionRelativeData.
      // It is here initialized in such a way that the sum of the n first positions is 1.0, where n is the
      // number of wall nodes.
      double w = 1.0 / mNeighbor_w->GetGeometry().size();
      array_1d<double, 4> weight;
      weight[0] = w;
      weight[1] = w;
      weight[2] = w;
      weight[3] = w;

      // Dummy variables: not used now
      double dummy1[3][3];
      DEM_SET_COMPONENTS_TO_ZERO_3x3(dummy1);
      array_1d<double, 3> dummy2 = ZeroVector(3);
      array_1d<double, 3> dummy3 = ZeroVector(3);
      int dummy4 = 0;

      mNeighbor_w->ComputeConditionRelativeData(mNeighborIndex, this, dummy1, distance, weight, dummy2, dummy3, dummy4);

      return distance;
    }
    else {
      return 0.0;
    }

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeDistanceToNeighborAdjusted() {
    KRATOS_TRY

    // Corrected distance based on a corrected contact radius

    if (mNeighborType & PARTICLE_NEIGHBOR) {
      const double r1 = GetParticleRadius();
      const double r2 = GetNeighborRadius();
      return sqrt(r1 * r1 - mContactRadiusAdjusted * mContactRadiusAdjusted) + sqrt(r2 * r2 - mContactRadiusAdjusted * mContactRadiusAdjusted);
    }
    else if (mNeighborType & WALL_NEIGHBOR) {
      const double r = GetParticleRadius();
      return sqrt(r * r - mContactRadiusAdjusted * mContactRadiusAdjusted);
    }
    else {
      return 0.0;
    }

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeSeparationToNeighbor() {
    KRATOS_TRY

    const double particle_radius = GetParticleRadius();
    const double neighbor_radius = GetNeighborRadius(); // must be zero for walls
    return mNeighborDistance - particle_radius - neighbor_radius;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeSeparationToNeighborAdjusted() {
    KRATOS_TRY

    const double particle_radius = GetParticleRadius();
    const double neighbor_radius = GetNeighborRadius(); // must be zero for walls
    return mNeighborDistanceAdjusted - particle_radius - neighbor_radius;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeFourierNumber() {
    KRATOS_TRY

    const double col_time_max = ComputeMaxCollisionTime();
    const double Rc_max       = ComputeMaxContactRadius();

    // Compute particle Fourier number
    double Fo_particle;
    if (Rc_max > 0.0)
      Fo_particle = GetParticleConductivity() * col_time_max / (GetParticleDensity() * GetParticleHeatCapacity() * Rc_max * Rc_max);
    else
      Fo_particle = 0.0;

    // Compute neighbor Fourier number
    double Fo_neighbor;

    if (mNeighborType & PARTICLE_NEIGHBOR) {
      if (Rc_max > 0.0)
        Fo_neighbor = GetNeighborConductivity() * col_time_max / (GetNeighborDensity() * GetNeighborHeatCapacity() * Rc_max * Rc_max);
      else
        Fo_neighbor = 0.0;
    }
    else {
      Fo_neighbor = Fo_particle;
    }

    // Assumption: average of both particles (only particle if neighbor is a wall)
    return (Fo_particle + Fo_neighbor) / 2.0;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeMaxCollisionTime() {
    KRATOS_TRY

    const double eff_radius             = ComputeEffectiveRadius();
    const double eff_mass               = ComputeEffectiveMass();
    const double eff_young              = ComputeEffectiveYoungReal(); // ATTENTION: Assumption: Original model was not assumed real Young modulus!
    const double impact_normal_velocity = fabs(GetContactParameters().impact_velocity[0]);

    if (impact_normal_velocity != 0.0)
      return 2.87 * pow(eff_mass * eff_mass / (eff_radius * eff_young * eff_young * impact_normal_velocity), 0.2);
    else
      return std::numeric_limits<double>::max();

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeMaxContactRadius() {
    KRATOS_TRY

    const double eff_radius             = ComputeEffectiveRadius();
    const double eff_mass               = ComputeEffectiveMass();
    const double eff_young              = ComputeEffectiveYoungReal(); // ATTENTION: Assumption: Original model was not assumed real Young modulus!
    const double impact_normal_velocity = fabs(GetContactParameters().impact_velocity[0]);

    return pow(15.0 * eff_mass * eff_radius * eff_radius * impact_normal_velocity * impact_normal_velocity / (16.0 * eff_young), 0.2);
    
    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeContactRadius() {
    KRATOS_TRY
    
    double Rc = 0.0;

    if (mNeighborInContact) {
      if (mNeighborType & PARTICLE_NEIGHBOR) {
        const double r1 = GetParticleRadius();
        const double r2 = GetNeighborRadius();
        Rc = sqrt(fabs(r1 * r1 - pow(((r1 * r1 - r2 * r2 + mNeighborDistance * mNeighborDistance) / (2.0 * mNeighborDistance)), 2.0)));
      }
      else if (mNeighborType & WALL_NEIGHBOR) {
        const double r = GetParticleRadius();
        Rc = sqrt(r * r - mNeighborDistance * mNeighborDistance);
      }
    }

    return Rc;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeEffectiveRadius() {
    KRATOS_TRY

    if (mNeighborType & PARTICLE_NEIGHBOR) {
      const double particle_radius = GetParticleRadius();
      const double neighbor_radius = GetNeighborRadius();
      return particle_radius * neighbor_radius / (particle_radius + neighbor_radius);
    }
    else if (mNeighborType & WALL_NEIGHBOR) {
      return GetParticleRadius();
    }
    else {
      return 0.0;
    }

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeEffectiveMass() {
    KRATOS_TRY

    if (mNeighborType & PARTICLE_NEIGHBOR) {
      const double particle_mass = GetParticleMass();
      const double neighbor_mass = GetNeighborMass();
      return particle_mass * neighbor_mass / (particle_mass + neighbor_mass);
    }
    else if (mNeighborType & WALL_NEIGHBOR) {
      return GetParticleMass();
    }
    else {
      return 0.0;
    }

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeEffectiveYoung() {
    KRATOS_TRY

    const double particle_young   = GetParticleYoung();
    const double particle_poisson = GetParticlePoisson();
    const double neighbor_young   = GetNeighborYoung();
    const double neighbor_poisson = GetNeighborPoisson();

    return 1.0 / ((1.0 - particle_poisson * particle_poisson) / particle_young + (1.0 - neighbor_poisson * neighbor_poisson) / neighbor_young);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeEffectiveYoungReal() {
    KRATOS_TRY

    const double particle_young   = GetParticleYoung() * mRealYoungRatio;
    const double particle_poisson = GetParticlePoisson();
    const double neighbor_young   = GetNeighborYoung() * mRealYoungRatio;
    const double neighbor_poisson = GetNeighborPoisson();

    return 1.0 / ((1.0 - particle_poisson * particle_poisson) / particle_young + (1.0 - neighbor_poisson * neighbor_poisson) / neighbor_young);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeEffectiveConductivity() {
    KRATOS_TRY

    const double particle_conductivity = GetParticleConductivity();
    const double neighbor_conductivity = GetNeighborConductivity();

    return particle_conductivity * neighbor_conductivity / (particle_conductivity + neighbor_conductivity);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeAverageConductivity() {
    KRATOS_TRY

    const double r1 = GetParticleRadius();
    const double r2 = GetNeighborRadius(); // must be zero for walls
    const double k1 = GetParticleConductivity();
    const double k2 = GetNeighborConductivity();

    return (r1 + r2) / (r1 / k1 + r2 / k2);

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Get/Set methods

  //------------------------------------------------------------------------------------------------------------
  ThermalDEMIntegrationScheme& ThermalSphericParticle::GetThermalIntegrationScheme() {
    return *mpThermalIntegrationScheme;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetYoung() {
    if (GetProperties().HasTable(TEMPERATURE, YOUNG_MODULUS)) {
      const auto& r_table = GetProperties().GetTable(TEMPERATURE, YOUNG_MODULUS);
      return r_table.GetValue(GetParticleTemperature());
    }
    else {
      return GetFastProperties()->GetYoung();
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetPoisson() {
    if (GetProperties().HasTable(TEMPERATURE, POISSON_RATIO)) {
      const auto& r_table = GetProperties().GetTable(TEMPERATURE, POISSON_RATIO);
      return r_table.GetValue(GetParticleTemperature());
    }
    else {
      return GetFastProperties()->GetPoisson();
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetDensity() {
    if (GetProperties().HasTable(TEMPERATURE, PARTICLE_DENSITY)) {
      const auto& r_table = GetProperties().GetTable(TEMPERATURE, PARTICLE_DENSITY);
      return r_table.GetValue(GetParticleTemperature());
    }
    else {
      return GetFastProperties()->GetDensity();
    }
  }

  //------------------------------------------------------------------------------------------------------------
  array_1d<double, 3> ThermalSphericParticle::GetParticleCoordinates() {
    return GetGeometry()[0].Coordinates();
  }

  //------------------------------------------------------------------------------------------------------------
  array_1d<double, 3> ThermalSphericParticle::GetParticleVelocity() {
    return GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleTemperature() {
    return GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleRadius() {
    return GetRadius();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleSurfaceArea() {
    return 4.0 * Globals::Pi * GetRadius() * GetRadius();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleCharacteristicLength() {
    return 2.0 * GetRadius();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleVolume() {
    return CalculateVolume();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleYoung() {
    return GetYoung();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticlePoisson() {
    return GetPoisson();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleDensity() {
    return GetDensity();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleMass() {
    return GetMass();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleHeatCapacity() {
    if (GetProperties().HasTable(TEMPERATURE, SPECIFIC_HEAT)) {
      const auto& r_table = GetProperties().GetTable(TEMPERATURE, SPECIFIC_HEAT);
      return r_table.GetValue(GetParticleTemperature());
    }
    else {
      return GetProperties()[SPECIFIC_HEAT]; // TODO: Use GetFastProperties?
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleConductivity() {
    if (GetProperties().HasTable(TEMPERATURE, THERMAL_CONDUCTIVITY)) {
      const auto& r_table = GetProperties().GetTable(TEMPERATURE, THERMAL_CONDUCTIVITY);
      return r_table.GetValue(GetParticleTemperature());
    }
    else {
      return GetProperties()[THERMAL_CONDUCTIVITY]; // TODO: Use GetFastProperties?
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleEmissivity() {
    if (GetProperties().HasTable(TEMPERATURE, EMISSIVITY)) {
      const auto& r_table = GetProperties().GetTable(TEMPERATURE, EMISSIVITY);
      return r_table.GetValue(GetParticleTemperature());
    }
    else {
      return GetProperties()[EMISSIVITY]; // TODO: Use GetFastProperties?
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleExpansionCoefficient() {
    if (GetProperties().HasTable(TEMPERATURE, THERMAL_EXPANSION_COEFFICIENT)) {
      const auto& r_table = GetProperties().GetTable(TEMPERATURE, THERMAL_EXPANSION_COEFFICIENT);
      return r_table.GetValue(GetParticleTemperature());
    }
    else {
      return GetProperties()[THERMAL_EXPANSION_COEFFICIENT]; // TODO: Use GetFastProperties?
    }
  }

  //------------------------------------------------------------------------------------------------------------
  array_1d<double, 3> ThermalSphericParticle::GetWallCoordinates() {
    return mNeighbor_w->GetGeometry().Center();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallTemperature() {
    // Assumption: wall temperature is the average of its nodes
    const double n_nodes = mNeighbor_w->GetGeometry().size();
    double wall_temp = 0.0;
    for (unsigned int i = 0; i < n_nodes; i++)
      wall_temp += mNeighbor_w->GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE);
    return wall_temp / n_nodes;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallRadius() {
    // Assumption: zero to be consistent with its use in formulations
    return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallYoung() {
    if (mNeighbor_w->GetProperties().HasTable(TEMPERATURE, YOUNG_MODULUS)) {
      const auto& r_table = mNeighbor_w->GetProperties().GetTable(TEMPERATURE, YOUNG_MODULUS);
      return r_table.GetValue(GetWallTemperature());
    }
    else {
      return mNeighbor_w->GetProperties()[YOUNG_MODULUS];
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallPoisson() {
    if (mNeighbor_w->GetProperties().HasTable(TEMPERATURE, POISSON_RATIO)) {
      const auto& r_table = mNeighbor_w->GetProperties().GetTable(TEMPERATURE, POISSON_RATIO);
      return r_table.GetValue(GetWallTemperature());
    }
    else {
      return mNeighbor_w->GetProperties()[POISSON_RATIO];
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallDensity() {
    if (mNeighbor_w->GetProperties().HasTable(TEMPERATURE, PARTICLE_DENSITY)) {
      const auto& r_table = mNeighbor_w->GetProperties().GetTable(TEMPERATURE, PARTICLE_DENSITY);
      return r_table.GetValue(GetWallTemperature());
    }
    else {
      return mNeighbor_w->GetProperties()[PARTICLE_DENSITY];
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallMass() {
    // Assumption: zero to be consistent with its use in formulations
    return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallHeatCapacity() {
    if (mNeighbor_w->GetProperties().HasTable(TEMPERATURE, SPECIFIC_HEAT)) {
      const auto& r_table = mNeighbor_w->GetProperties().GetTable(TEMPERATURE, SPECIFIC_HEAT);
      return r_table.GetValue(GetWallTemperature());
    }
    else {
      return mNeighbor_w->GetProperties()[SPECIFIC_HEAT];
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallConductivity() {
    if (mNeighbor_w->GetProperties().HasTable(TEMPERATURE, THERMAL_CONDUCTIVITY)) {
      const auto& r_table = mNeighbor_w->GetProperties().GetTable(TEMPERATURE, THERMAL_CONDUCTIVITY);
      return r_table.GetValue(GetWallTemperature());
    }
    else {
      return mNeighbor_w->GetProperties()[THERMAL_CONDUCTIVITY];
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallEmissivity() {
    if (mNeighbor_w->GetProperties().HasTable(TEMPERATURE, EMISSIVITY)) {
      const auto& r_table = mNeighbor_w->GetProperties().GetTable(TEMPERATURE, EMISSIVITY);
      return r_table.GetValue(GetWallTemperature());
    }
    else {
      return mNeighbor_w->GetProperties()[EMISSIVITY];
    }
  }

  //------------------------------------------------------------------------------------------------------------
  array_1d<double, 3> ThermalSphericParticle::GetNeighborCoordinates() {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleCoordinates();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallCoordinates();
    else
      return vector<double>();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborTemperature() {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleTemperature();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallTemperature();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborRadius() {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleRadius();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallRadius();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborYoung() {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleYoung();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallYoung();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborPoisson() {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticlePoisson();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallPoisson();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborDensity() {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleDensity();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallDensity();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborMass() {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleMass();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallMass();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborHeatCapacity() {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleHeatCapacity();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallHeatCapacity();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborConductivity() {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleConductivity();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallConductivity();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborEmissivity() {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleEmissivity();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallEmissivity();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetContactDynamicFrictionCoefficient() {
    if (mNeighborType & PARTICLE_NEIGHBOR) {
      Properties& properties_of_contact = GetProperties().GetSubProperties(mNeighbor_p->GetProperties().Id());
      return properties_of_contact[DYNAMIC_FRICTION];
    }
    else if (mNeighborType & WALL_NEIGHBOR) {
      Properties& properties_of_contact = GetProperties().GetSubProperties(mNeighbor_w->GetProperties().Id());
      return properties_of_contact[DYNAMIC_FRICTION];
    }
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  typename ThermalSphericParticle::ContactParams ThermalSphericParticle::GetContactParameters() {
    if (mNeighborType & PARTICLE_NEIGHBOR && mContactParamsParticle.count(mNeighbor_p)) {
      return mContactParamsParticle[mNeighbor_p];
    }
    else if (mNeighborType & WALL_NEIGHBOR && mContactParamsWall.count(mNeighbor_w)) {
      return mContactParamsWall[mNeighbor_w];
    }
    else {
      ContactParams null_param;
      null_param.updated_step = 0;
      null_param.impact_time  = 0.0;
      null_param.impact_velocity.assign(2, 0.0);
      null_param.local_velocity.assign(2, 0.0);
      null_param.local_force.assign(2, 0.0);
      return null_param;
    }
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetThermalIntegrationScheme(ThermalDEMIntegrationScheme::Pointer& scheme) {
    mpThermalIntegrationScheme = scheme->CloneRaw();
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetParticleTemperature(const double temperature) {
    GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE) = temperature;
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetParticleHeatFlux(const double heat_flux) {
    GetGeometry()[0].FastGetSolutionStepValue(HEATFLUX) = heat_flux;
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetParticlePrescribedHeatFluxSurface(const double heat_flux) {
    mPrescribedHeatFluxSurface = heat_flux;
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetParticlePrescribedHeatFluxVolume(const double heat_flux) {
    mPrescribedHeatFluxVolume = heat_flux;
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetParticleRadius(const double radius) {
    SetRadius(radius);
    GetGeometry()[0].FastGetSolutionStepValue(RADIUS) = radius;
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetParticleMass(const double mass) {
    SetMass(mass);
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetParticleMomentInertia(const double moment_inertia) {
    GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = moment_inertia;
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetParticleRealYoungRatio(const double ratio) {
    mRealYoungRatio = ratio;
  }

} // namespace Kratos
