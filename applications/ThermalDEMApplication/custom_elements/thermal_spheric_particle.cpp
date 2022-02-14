//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \
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
    mpThermalIntegrationScheme   = NULL;
    mpNumericalIntegrationMethod = NULL;
    mpDirectConductionModel      = NULL;
    mpIndirectConductionModel    = NULL;
    mpConvectionModel            = NULL;
    mpRadiationModel             = NULL;
    mpFrictionModel              = NULL;
  }

  ThermalSphericParticle::ThermalSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry):SphericParticle(NewId, pGeometry) {
    mpThermalIntegrationScheme   = NULL;
    mpNumericalIntegrationMethod = NULL;
    mpDirectConductionModel      = NULL;
    mpIndirectConductionModel    = NULL;
    mpConvectionModel            = NULL;
    mpRadiationModel             = NULL;
    mpFrictionModel              = NULL;
  }

  ThermalSphericParticle::ThermalSphericParticle(IndexType NewId, NodesArrayType const& ThisNodes):SphericParticle(NewId, ThisNodes) {
    mpThermalIntegrationScheme   = NULL;
    mpNumericalIntegrationMethod = NULL;
    mpDirectConductionModel      = NULL;
    mpIndirectConductionModel    = NULL;
    mpConvectionModel            = NULL;
    mpRadiationModel             = NULL;
    mpFrictionModel              = NULL;
  }

  ThermalSphericParticle::ThermalSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):SphericParticle(NewId, pGeometry, pProperties) {
    mpThermalIntegrationScheme   = NULL;
    mpNumericalIntegrationMethod = NULL;
    mpDirectConductionModel      = NULL;
    mpIndirectConductionModel    = NULL;
    mpConvectionModel            = NULL;
    mpRadiationModel             = NULL;
    mpFrictionModel              = NULL;
  }

  Element::Pointer ThermalSphericParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
    return SphericParticle::Pointer(new ThermalSphericParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
  }

  ThermalSphericParticle::~ThermalSphericParticle() {
    if (mpThermalIntegrationScheme != NULL) {
      delete mpThermalIntegrationScheme;
      mpThermalIntegrationScheme = NULL;
    }
    if (mpNumericalIntegrationMethod != NULL) {
      delete mpNumericalIntegrationMethod;
      mpNumericalIntegrationMethod = NULL;
    }
    if (mpDirectConductionModel != NULL) {
      delete mpDirectConductionModel;
      mpDirectConductionModel = NULL;
    }
    if (mpIndirectConductionModel != NULL) {
      delete mpIndirectConductionModel;
      mpIndirectConductionModel = NULL;
    }
    if (mpConvectionModel != NULL) {
      delete mpConvectionModel;
      mpConvectionModel = NULL;
    }
    if (mpRadiationModel != NULL) {
      delete mpRadiationModel;
      mpRadiationModel = NULL;
    }
    if (mpFrictionModel != NULL) {
      delete mpFrictionModel;
      mpFrictionModel = NULL;
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

    // Set time integration scheme
    ThermalDEMIntegrationScheme::Pointer& thermal_integration_scheme = GetProperties()[THERMAL_INTEGRATION_SCHEME_POINTER];
    SetThermalIntegrationScheme(thermal_integration_scheme);

    // Set numerical integration method
    NumericalIntegrationMethod::Pointer& numerical_integration_method = GetProperties()[NUMERICAL_INTEGRATION_METHOD_POINTER];
    SetNumericalIntegrationMethod(numerical_integration_method);

    // Set constitutive laws (thermal models)
    HeatExchangeMechanism::Pointer& direct_conduction_model = GetProperties()[DIRECT_CONDUCTION_MODEL_POINTER];
    SetDirectConductionModel(direct_conduction_model);

    HeatExchangeMechanism::Pointer& indirect_conduction_model = GetProperties()[INDIRECT_CONDUCTION_MODEL_POINTER];
    SetIndirectConductionModel(indirect_conduction_model);

    HeatExchangeMechanism::Pointer& convection_model = GetProperties()[CONVECTION_MODEL_POINTER];
    SetConvectionModel(convection_model);

    HeatExchangeMechanism::Pointer& radiation_model = GetProperties()[RADIATION_MODEL_POINTER];
    SetRadiationModel(radiation_model);

    HeatGenerationMechanism::Pointer& friction_model = GetProperties()[FRICTION_MODEL_POINTER];
    SetFrictionModel(friction_model);

    // Set flag to store contact parameters during mechanical loop over neighbors
    mStoreContactParam = this->Is(DEMThermalFlags::HAS_MOTION)        &&
                        (this->Is(DEMThermalFlags::HAS_FRICTION_HEAT) ||
                        (this->Is(DEMThermalFlags::HAS_DIRECT_CONDUCTION) && r_process_info[DIRECT_CONDUCTION_MODEL_NAME].compare("collisional") == 0));    

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

      // ATTENTION:
      // Maximum of 2 contact neighbor walls (maybe can be removed).
      // A rare bug was observed in the particle-wall contact in which 3 neighbor walls were detected when there was only 1.
      // To avoid overcomputing heat transfer, this limit of 2 neighbors is being imposed.
      if (i > 1)
        return;
    }

    // Compute heat fluxes with noncontact neighbor walls
    // ATTENTION:
    // Only one element of a noncontact wall is elected to represent the entire wall,
    // which is the edge(2D)/face(3D) intersected by the normal vector between the wall and the particle.
    // This elected element represents an infinity wall.
    for (unsigned int i = 0; i < mNeighbourNonContactRigidFaces.size(); i++) {
      if (mNeighbourNonContactRigidFaces[i] == NULL) continue;
      mNeighbor_w    = dynamic_cast<DEMWall*>(mNeighbourNonContactRigidFaces[i]);
      mNeighborType  = WALL_NEIGHBOR_NONCONTACT;
      mNeighborIndex = i;
      ComputeHeatFluxWithNeighbor(r_process_info);
    }

    // Finalize radiation computation of continuous methods
    if (this->Is(DEMThermalFlags::HAS_RADIATION))
      mRadiationHeatFlux += GetRadiationModel().FinalizeHeatFlux(r_process_info, this);

    // Compute convection with surrounding fluid
    if (this->Is(DEMThermalFlags::HAS_CONVECTION))
      mConvectionHeatFlux += GetConvectionModel().ComputeHeatFlux(r_process_info, this);

    // Prescribed heat flux over surface area
    if (mPrescribedHeatFluxSurface != 0.0)
      mPrescribedHeatFlux += mPrescribedHeatFluxSurface * GetParticleSurfaceArea();

    // Prescribed heat source over volume
    if (mPrescribedHeatFluxVolume != 0.0)
      mPrescribedHeatFlux += mPrescribedHeatFluxVolume * GetParticleVolume();

    // Sum up heat fluxes contributions
    mTotalHeatFlux = mConductionDirectHeatFlux + mConductionIndirectHeatFlux + mRadiationHeatFlux + mFrictionHeatFlux + mConvectionHeatFlux + mPrescribedHeatFlux;
    SetParticleHeatFlux(mTotalHeatFlux);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputeHeatFluxWithNeighbor(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Check if neighbor is adiabatic
    if (CheckAdiabaticNeighbor())
      return;

    // Compute simulated or adjusted interaction properties
    ComputeInteractionProps(r_process_info);

    // Heat transfer mechanisms
    if (this->Is(DEMThermalFlags::HAS_DIRECT_CONDUCTION))
      mConductionDirectHeatFlux += GetDirectConductionModel().ComputeHeatFlux(r_process_info, this);

    if (this->Is(DEMThermalFlags::HAS_INDIRECT_CONDUCTION))
      mConductionIndirectHeatFlux += GetIndirectConductionModel().ComputeHeatFlux(r_process_info, this);

    if (this->Is(DEMThermalFlags::HAS_RADIATION))
      mRadiationHeatFlux += GetRadiationModel().ComputeHeatFlux(r_process_info, this);

    if (this->Is(DEMThermalFlags::HAS_FRICTION_HEAT) && this->Is(DEMThermalFlags::HAS_MOTION))
      mFrictionHeatFlux += GetFrictionModel().ComputeHeatGeneration(r_process_info, this);

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
      std::string model = r_process_info[ADJUSTED_CONTACT_MODEL_NAME];
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

    // Update temperature dependent properties
    if (mIsTimeToSolve)
      UpdateTemperatureDependentRadius(r_process_info);

    KRATOS_CATCH("")
  }

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
  // Contact adjustment models

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
    const double identation     = std::max(-mNeighborSeparation, 0.0);

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
  // Integration expressions

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::EvalIntegrandSurrLayer(NumericalIntegrationMethod* method) {
    KRATOS_TRY

    const double r    = method->mCoord;
    const double d    = method->mParams.p1;
    const double dmin = method->mParams.p2;
    const double r1   = method->mParams.p3;
    const double r2   = method->mParams.p4;

    return 2.0 * Globals::Pi * r / std::max(dmin, d - sqrt(r1 * r1 - r * r) - sqrt(r2 * r2 - r * r));

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::EvalIntegrandVoronoiWall(NumericalIntegrationMethod* method) {
    KRATOS_TRY

    const double r   = method->mCoord;
    const double d   = method->mParams.p1;
    const double kf  = method->mParams.p2;
    const double kp  = method->mParams.p3;
    const double rp  = method->mParams.p4;
    const double rij = method->mParams.p5;

    return 2.0 * Globals::Pi * r / ((sqrt(rp * rp - r * r) - r * d / rij) / kp + (d - sqrt(rp * rp - r * r)) / kf);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::EvalIntegrandVoronoiMono(NumericalIntegrationMethod* method) {
    KRATOS_TRY

    const double r    = method->mCoord;
    const double d    = method->mParams.p1;
    const double kf   = method->mParams.p2;
    const double keff = method->mParams.p3;
    const double rp   = method->mParams.p4;
    const double rij  = method->mParams.p5;

    return 2.0 * Globals::Pi * r / ((sqrt(rp * rp - r * r) - r * d / (2.0 * rij)) / keff + 2.0 * (d / 2.0 - sqrt(rp * rp - r * r)) / kf);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::EvalIntegrandVoronoiMulti(NumericalIntegrationMethod* method) {
    KRATOS_TRY

    const double r    = method->mCoord;
    const double d    = method->mParams.p1;
    const double kf   = method->mParams.p2;
    const double k1   = method->mParams.p3;
    const double k2   = method->mParams.p4;
    const double r1   = method->mParams.p5;
    const double r2   = method->mParams.p6;
    const double rij  = method->mParams.p7;
    const double rij_ = method->mParams.p8;
    const double D1   = method->mParams.p9;
    const double D2   = method->mParams.p10;

    const double beta1 = sqrt(r1 * r1 - r * r);
    const double beta2 = sqrt(r2 * r2 - r * r);

    return 2.0 * Globals::Pi * r / ((beta1 - D1 * r / rij) / k1 + (beta2 - D2 * r / rij_) / k2 + (d - beta1 - beta2) / kf);

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Auxiliary computations

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputeAddedSearchDistance(const ProcessInfo& r_process_info, double& added_search_distance) {
    KRATOS_TRY

    if (this->Is(DEMThermalFlags::HAS_INDIRECT_CONDUCTION)) {
      const double model_search_distance  = GetIndirectConductionModel().GetSearchDistance(r_process_info, this);
      const double current_added_distance = added_search_distance;
      added_search_distance = std::max(current_added_distance, model_search_distance);
    }

    if (this->Is(DEMThermalFlags::HAS_RADIATION)) {
      const double model_search_distance  = GetRadiationModel().GetSearchDistance(r_process_info, this);
      const double current_added_distance = added_search_distance;
      added_search_distance = std::max(current_added_distance, model_search_distance);
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
    if (r_process_info[VORONOI_METHOD_NAME].compare("tesselation") == 0) {
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
    else if (r_process_info[VORONOI_METHOD_NAME].compare("porosity") == 0) {
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
  bool ThermalSphericParticle::CheckAdiabaticNeighbor(void) {
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
  double ThermalSphericParticle::ComputeDistanceToNeighbor(void) {
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
  double ThermalSphericParticle::ComputeDistanceToNeighborAdjusted(void) {
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
  double ThermalSphericParticle::ComputeSeparationToNeighbor(void) {
    KRATOS_TRY

    const double particle_radius = GetParticleRadius();
    const double neighbor_radius = GetNeighborRadius(); // must be zero for walls
    return mNeighborDistance - particle_radius - neighbor_radius;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeSeparationToNeighborAdjusted(void) {
    KRATOS_TRY

    const double particle_radius = GetParticleRadius();
    const double neighbor_radius = GetNeighborRadius(); // must be zero for walls
    return mNeighborDistanceAdjusted - particle_radius - neighbor_radius;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeFourierNumber(void) {
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
  double ThermalSphericParticle::ComputeMaxCollisionTime(void) {
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
  double ThermalSphericParticle::ComputeMaxContactRadius(void) {
    KRATOS_TRY

    const double eff_radius             = ComputeEffectiveRadius();
    const double eff_mass               = ComputeEffectiveMass();
    const double eff_young              = ComputeEffectiveYoungReal(); // ATTENTION: Assumption: Original model was not assumed real Young modulus!
    const double impact_normal_velocity = fabs(GetContactParameters().impact_velocity[0]);

    return pow(15.0 * eff_mass * eff_radius * eff_radius * impact_normal_velocity * impact_normal_velocity / (16.0 * eff_young), 0.2);
    
    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeContactRadius(void) {
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
  double ThermalSphericParticle::ComputeEffectiveRadius(void) {
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
  double ThermalSphericParticle::ComputeEffectiveMass(void) {
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
  double ThermalSphericParticle::ComputeEffectiveYoung(void) {
    KRATOS_TRY

    const double particle_young   = GetParticleYoung();
    const double particle_poisson = GetParticlePoisson();
    const double neighbor_young   = GetNeighborYoung();
    const double neighbor_poisson = GetNeighborPoisson();

    return 1.0 / ((1.0 - particle_poisson * particle_poisson) / particle_young + (1.0 - neighbor_poisson * neighbor_poisson) / neighbor_young);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeEffectiveYoungReal(void) {
    KRATOS_TRY

    const double particle_young   = GetParticleYoung() * mRealYoungRatio;
    const double particle_poisson = GetParticlePoisson();
    const double neighbor_young   = GetNeighborYoung() * mRealYoungRatio;
    const double neighbor_poisson = GetNeighborPoisson();

    return 1.0 / ((1.0 - particle_poisson * particle_poisson) / particle_young + (1.0 - neighbor_poisson * neighbor_poisson) / neighbor_young);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeEffectiveConductivity(void) {
    KRATOS_TRY

    const double particle_conductivity = GetParticleConductivity();
    const double neighbor_conductivity = GetNeighborConductivity();

    return particle_conductivity * neighbor_conductivity / (particle_conductivity + neighbor_conductivity);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeAverageConductivity(void) {
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
  ThermalDEMIntegrationScheme& ThermalSphericParticle::GetThermalIntegrationScheme(void) {
    return *mpThermalIntegrationScheme;
  }

  //------------------------------------------------------------------------------------------------------------
  NumericalIntegrationMethod& ThermalSphericParticle::GetNumericalIntegrationMethod(void) {
    return *mpNumericalIntegrationMethod;
  }

  //------------------------------------------------------------------------------------------------------------
  HeatExchangeMechanism& ThermalSphericParticle::GetDirectConductionModel(void) {
    return *mpDirectConductionModel;
  }

  //------------------------------------------------------------------------------------------------------------
  HeatExchangeMechanism& ThermalSphericParticle::GetIndirectConductionModel(void) {
    return *mpIndirectConductionModel;
  }

  //------------------------------------------------------------------------------------------------------------
  HeatExchangeMechanism& ThermalSphericParticle::GetConvectionModel(void) {
    return *mpConvectionModel;
  }

  //------------------------------------------------------------------------------------------------------------
  HeatExchangeMechanism& ThermalSphericParticle::GetRadiationModel(void) {
    return *mpRadiationModel;
  }

  //------------------------------------------------------------------------------------------------------------
  HeatGenerationMechanism& ThermalSphericParticle::GetFrictionModel(void) {
    return *mpFrictionModel;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetYoung(void) {
    if (GetProperties().HasTable(TEMPERATURE, YOUNG_MODULUS)) {
      const auto& r_table = GetProperties().GetTable(TEMPERATURE, YOUNG_MODULUS);
      return r_table.GetValue(GetParticleTemperature());
    }
    else {
      return GetFastProperties()->GetYoung();
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetPoisson(void) {
    if (GetProperties().HasTable(TEMPERATURE, POISSON_RATIO)) {
      const auto& r_table = GetProperties().GetTable(TEMPERATURE, POISSON_RATIO);
      return r_table.GetValue(GetParticleTemperature());
    }
    else {
      return GetFastProperties()->GetPoisson();
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetDensity(void) {
    if (GetProperties().HasTable(TEMPERATURE, PARTICLE_DENSITY)) {
      const auto& r_table = GetProperties().GetTable(TEMPERATURE, PARTICLE_DENSITY);
      return r_table.GetValue(GetParticleTemperature());
    }
    else {
      return GetFastProperties()->GetDensity();
    }
  }

  //------------------------------------------------------------------------------------------------------------
  array_1d<double, 3> ThermalSphericParticle::GetParticleCoordinates(void) {
    return GetGeometry()[0].Coordinates();
  }

  //------------------------------------------------------------------------------------------------------------
  array_1d<double, 3> ThermalSphericParticle::GetParticleVelocity(void) {
    return GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleTemperature(void) {
    return GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleRadius(void) {
    return GetRadius();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleSurfaceArea(void) {
    return 4.0 * Globals::Pi * GetRadius() * GetRadius();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleCharacteristicLength(void) {
    return 2.0 * GetRadius();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleVolume(void) {
    return CalculateVolume();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleYoung(void) {
    return GetYoung();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticlePoisson(void) {
    return GetPoisson();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleDensity(void) {
    return GetDensity();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleMass(void) {
    return GetMass();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleHeatCapacity(void) {
    if (GetProperties().HasTable(TEMPERATURE, SPECIFIC_HEAT)) {
      const auto& r_table = GetProperties().GetTable(TEMPERATURE, SPECIFIC_HEAT);
      return r_table.GetValue(GetParticleTemperature());
    }
    else {
      return GetProperties()[SPECIFIC_HEAT]; // TODO: Use GetFastProperties?
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleConductivity(void) {
    if (GetProperties().HasTable(TEMPERATURE, THERMAL_CONDUCTIVITY)) {
      const auto& r_table = GetProperties().GetTable(TEMPERATURE, THERMAL_CONDUCTIVITY);
      return r_table.GetValue(GetParticleTemperature());
    }
    else {
      return GetProperties()[THERMAL_CONDUCTIVITY]; // TODO: Use GetFastProperties?
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleEmissivity(void) {
    if (GetProperties().HasTable(TEMPERATURE, EMISSIVITY)) {
      const auto& r_table = GetProperties().GetTable(TEMPERATURE, EMISSIVITY);
      return r_table.GetValue(GetParticleTemperature());
    }
    else {
      return GetProperties()[EMISSIVITY]; // TODO: Use GetFastProperties?
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleExpansionCoefficient(void) {
    if (GetProperties().HasTable(TEMPERATURE, THERMAL_EXPANSION_COEFFICIENT)) {
      const auto& r_table = GetProperties().GetTable(TEMPERATURE, THERMAL_EXPANSION_COEFFICIENT);
      return r_table.GetValue(GetParticleTemperature());
    }
    else {
      return GetProperties()[THERMAL_EXPANSION_COEFFICIENT]; // TODO: Use GetFastProperties?
    }
  }

  //------------------------------------------------------------------------------------------------------------
  array_1d<double, 3> ThermalSphericParticle::GetWallCoordinates(void) {
    return mNeighbor_w->GetGeometry().Center();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallTemperature(void) {
    // Assumption: wall temperature is the average of its nodes
    const double n_nodes = mNeighbor_w->GetGeometry().size();
    double wall_temp = 0.0;
    for (unsigned int i = 0; i < n_nodes; i++)
      wall_temp += mNeighbor_w->GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE);
    return wall_temp / n_nodes;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallRadius(void) {
    // Assumption: zero to be consistent with its use in formulations
    return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallSurfaceArea(void) {
    // Assumption: zero
    return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallYoung(void) {
    if (mNeighbor_w->GetProperties().HasTable(TEMPERATURE, YOUNG_MODULUS)) {
      const auto& r_table = mNeighbor_w->GetProperties().GetTable(TEMPERATURE, YOUNG_MODULUS);
      return r_table.GetValue(GetWallTemperature());
    }
    else {
      return mNeighbor_w->GetProperties()[YOUNG_MODULUS];
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallPoisson(void) {
    if (mNeighbor_w->GetProperties().HasTable(TEMPERATURE, POISSON_RATIO)) {
      const auto& r_table = mNeighbor_w->GetProperties().GetTable(TEMPERATURE, POISSON_RATIO);
      return r_table.GetValue(GetWallTemperature());
    }
    else {
      return mNeighbor_w->GetProperties()[POISSON_RATIO];
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallDensity(void) {
    if (mNeighbor_w->GetProperties().HasTable(TEMPERATURE, PARTICLE_DENSITY)) {
      const auto& r_table = mNeighbor_w->GetProperties().GetTable(TEMPERATURE, PARTICLE_DENSITY);
      return r_table.GetValue(GetWallTemperature());
    }
    else {
      return mNeighbor_w->GetProperties()[PARTICLE_DENSITY];
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallMass(void) {
    // Assumption: zero to be consistent with its use in formulations
    return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallHeatCapacity(void) {
    if (mNeighbor_w->GetProperties().HasTable(TEMPERATURE, SPECIFIC_HEAT)) {
      const auto& r_table = mNeighbor_w->GetProperties().GetTable(TEMPERATURE, SPECIFIC_HEAT);
      return r_table.GetValue(GetWallTemperature());
    }
    else {
      return mNeighbor_w->GetProperties()[SPECIFIC_HEAT];
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallConductivity(void) {
    if (mNeighbor_w->GetProperties().HasTable(TEMPERATURE, THERMAL_CONDUCTIVITY)) {
      const auto& r_table = mNeighbor_w->GetProperties().GetTable(TEMPERATURE, THERMAL_CONDUCTIVITY);
      return r_table.GetValue(GetWallTemperature());
    }
    else {
      return mNeighbor_w->GetProperties()[THERMAL_CONDUCTIVITY];
    }
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetWallEmissivity(void) {
    if (mNeighbor_w->GetProperties().HasTable(TEMPERATURE, EMISSIVITY)) {
      const auto& r_table = mNeighbor_w->GetProperties().GetTable(TEMPERATURE, EMISSIVITY);
      return r_table.GetValue(GetWallTemperature());
    }
    else {
      return mNeighbor_w->GetProperties()[EMISSIVITY];
    }
  }

  //------------------------------------------------------------------------------------------------------------
  array_1d<double, 3> ThermalSphericParticle::GetNeighborCoordinates(void) {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleCoordinates();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallCoordinates();
    else
      return vector<double>();
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborTemperature(void) {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleTemperature();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallTemperature();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborRadius(void) {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleRadius();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallRadius();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborSurfaceArea(void) {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleSurfaceArea();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallSurfaceArea();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborYoung(void) {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleYoung();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallYoung();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborPoisson(void) {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticlePoisson();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallPoisson();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborDensity(void) {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleDensity();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallDensity();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborMass(void) {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleMass();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallMass();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborHeatCapacity(void) {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleHeatCapacity();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallHeatCapacity();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborConductivity(void) {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleConductivity();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallConductivity();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetNeighborEmissivity(void) {
    if (mNeighborType & PARTICLE_NEIGHBOR)
      return mNeighbor_p->GetParticleEmissivity();
    else if (mNeighborType & WALL_NEIGHBOR)
      return GetWallEmissivity();
    else
      return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetContactDynamicFrictionCoefficient(void) {
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
  typename ThermalSphericParticle::ContactParams ThermalSphericParticle::GetContactParameters(void) {
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
  void ThermalSphericParticle::SetNumericalIntegrationMethod(NumericalIntegrationMethod::Pointer& method) {
    mpNumericalIntegrationMethod = method->CloneRaw();
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetDirectConductionModel(HeatExchangeMechanism::Pointer& model) {
    mpDirectConductionModel = model->CloneRaw();
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetIndirectConductionModel(HeatExchangeMechanism::Pointer& model) {
    mpIndirectConductionModel = model->CloneRaw();
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetConvectionModel(HeatExchangeMechanism::Pointer& model) {
    mpConvectionModel = model->CloneRaw();
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetRadiationModel(HeatExchangeMechanism::Pointer& model) {
    mpRadiationModel = model->CloneRaw();
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetFrictionModel(HeatGenerationMechanism::Pointer& model) {
    mpFrictionModel = model->CloneRaw();
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
