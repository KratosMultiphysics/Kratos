//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
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
    mpGenerationModel            = NULL;
    mpRealContactModel           = NULL;
  }

  ThermalSphericParticle::ThermalSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry):SphericParticle(NewId, pGeometry) {
    mpThermalIntegrationScheme   = NULL;
    mpNumericalIntegrationMethod = NULL;
    mpDirectConductionModel      = NULL;
    mpIndirectConductionModel    = NULL;
    mpConvectionModel            = NULL;
    mpRadiationModel             = NULL;
    mpGenerationModel            = NULL;
    mpRealContactModel           = NULL;
  }

  ThermalSphericParticle::ThermalSphericParticle(IndexType NewId, NodesArrayType const& ThisNodes):SphericParticle(NewId, ThisNodes) {
    mpThermalIntegrationScheme   = NULL;
    mpNumericalIntegrationMethod = NULL;
    mpDirectConductionModel      = NULL;
    mpIndirectConductionModel    = NULL;
    mpConvectionModel            = NULL;
    mpRadiationModel             = NULL;
    mpGenerationModel            = NULL;
    mpRealContactModel           = NULL;
  }

  ThermalSphericParticle::ThermalSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):SphericParticle(NewId, pGeometry, pProperties) {
    mpThermalIntegrationScheme   = NULL;
    mpNumericalIntegrationMethod = NULL;
    mpDirectConductionModel      = NULL;
    mpIndirectConductionModel    = NULL;
    mpConvectionModel            = NULL;
    mpRadiationModel             = NULL;
    mpGenerationModel            = NULL;
    mpRealContactModel           = NULL;
  }

  Element::Pointer ThermalSphericParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
    return Element::Pointer(new ThermalSphericParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
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
    if (mpGenerationModel != NULL) {
      delete mpGenerationModel;
      mpGenerationModel = NULL;
    }
    if (mpRealContactModel != NULL) {
      delete mpRealContactModel;
      mpRealContactModel = NULL;
    }
  }

  //=====================================================================================================================================================================================
  // Initialization methods

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::Initialize(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    Properties& r_properties = GetProperties();

    // Dimension
    mDimension = r_process_info[DOMAIN_SIZE];

    // Initialize base class
    SphericParticle::Initialize(r_process_info);

    // Set pointers to to auxiliary objects
    HeatExchangeMechanism::Pointer&       direct_conduction_model      = r_properties[DIRECT_CONDUCTION_MODEL_POINTER];
    HeatExchangeMechanism::Pointer&       indirect_conduction_model    = r_properties[INDIRECT_CONDUCTION_MODEL_POINTER];
    HeatExchangeMechanism::Pointer&       convection_model             = r_properties[CONVECTION_MODEL_POINTER];
    HeatExchangeMechanism::Pointer&       radiation_model              = r_properties[RADIATION_MODEL_POINTER];
    HeatGenerationMechanism::Pointer&     generation_model             = r_properties[GENERATION_MODEL_POINTER];
    RealContactModel::Pointer&            real_contact_model           = r_properties[REAL_CONTACT_MODEL_POINTER];
    ThermalDEMIntegrationScheme::Pointer& thermal_integration_scheme   = r_properties[THERMAL_INTEGRATION_SCHEME_POINTER];
    NumericalIntegrationMethod::Pointer&  numerical_integration_method = r_properties[NUMERICAL_INTEGRATION_METHOD_POINTER];

    SetDirectConductionModel(direct_conduction_model);
    SetIndirectConductionModel(indirect_conduction_model);
    SetConvectionModel(convection_model);
    SetRadiationModel(radiation_model);
    SetGenerationModel(generation_model);
    SetRealContactModel(real_contact_model);
    SetThermalIntegrationScheme(thermal_integration_scheme);
    SetNumericalIntegrationMethod(numerical_integration_method);

    // Set flags
    mHasMotion = r_process_info[MOTION_OPTION];
    
    mHasVariableRadius = (r_properties.Has(THERMAL_EXPANSION_COEFFICIENT) && r_properties[THERMAL_EXPANSION_COEFFICIENT] != 0.0) ||
                          r_properties.HasTable(TEMPERATURE, THERMAL_EXPANSION_COEFFICIENT);

    mStoreContactParam = mHasMotion &&
                        (r_process_info[HEAT_GENERATION_OPTION]  ||
                        (r_process_info[DIRECT_CONDUCTION_OPTION] && r_process_info[DIRECT_CONDUCTION_MODEL_NAME].compare("collisional") == 0) || 
                        (r_process_info[REAL_CONTACT_OPTION]));
    
    // Clear maps
    mContactParamsParticle.clear();
    mContactParamsWall.clear();

    // Initialize accumulated energy dissipations
    mPreviousViscodampingEnergy            = 0.0;
    mPreviousFrictionalEnergy              = 0.0;
    mPreviousRollResistEnergy              = 0.0;
    mGenerationThermalEnergy_damp_particle = 0.0;
    mGenerationThermalEnergy_damp_wall     = 0.0;
    mGenerationThermalEnergy_slid_particle = 0.0;
    mGenerationThermalEnergy_slid_wall     = 0.0;
    mGenerationThermalEnergy_roll_particle = 0.0;
    mGenerationThermalEnergy_roll_wall     = 0.0;
    mGenerationThermalEnergy_wall_annul    = 0.0;
    mGenerationThermalEnergy_wall_sides    = 0.0;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::InitializeSolutionStep(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Initialize base class
    if (mHasMotion)
      SphericParticle::InitializeSolutionStep(r_process_info);

    // Check if it is time to evaluate thermal problem
    const int step = r_process_info[TIME_STEPS];
    const int freq = r_process_info[THERMAL_FREQUENCY];
    mIsTimeToSolve = (step > 0) && (freq != 0) && (step - 1) % freq == 0;

    // Number of steps passed between thermal evaluation steps
    mNumStepsEval = (step == 1) ? 1 : freq;

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
    mConductionDirectHeatFlux         = 0.0;
    mConductionIndirectHeatFlux       = 0.0;
    mConvectionHeatFlux               = 0.0;
    mRadiationHeatFlux                = 0.0;
    mGenerationHeatFlux               = 0.0;
    mGenerationHeatFlux_damp_particle = 0.0;
    mGenerationHeatFlux_damp_wall     = 0.0;
    mGenerationHeatFlux_slid_particle = 0.0;
    mGenerationHeatFlux_slid_wall     = 0.0;
    mGenerationHeatFlux_roll_particle = 0.0;
    mGenerationHeatFlux_roll_wall     = 0.0;
    mPrescribedHeatFlux               = 0.0;
    mTotalHeatFlux                    = 0.0;

    // Initialize environment-related variables for radiation
    if (r_process_info[RADIATION_OPTION]) {
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
    if (mHasMotion)
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
    if (r_process_info[RADIATION_OPTION])
      mRadiationHeatFlux += GetRadiationModel().FinalizeHeatFlux(r_process_info, this);

    // Convection with surrounding fluid
    if (r_process_info[CONVECTION_OPTION])
      mConvectionHeatFlux += GetConvectionModel().ComputeHeatFlux(r_process_info, this);

    // Prescribed heat flux over surface area
    if (mPrescribedHeatFluxSurface != 0.0)
      mPrescribedHeatFlux += mPrescribedHeatFluxSurface * GetParticleSurfaceArea();

    // Prescribed heat source over volume
    if (mPrescribedHeatFluxVolume != 0.0)
      mPrescribedHeatFlux += mPrescribedHeatFluxVolume * GetParticleVolume();

    // Sum up heat fluxes contributions
    mTotalHeatFlux = mConductionDirectHeatFlux + mConductionIndirectHeatFlux + mRadiationHeatFlux + mConvectionHeatFlux + mGenerationHeatFlux + mPrescribedHeatFlux;
    SetParticleHeatFlux(mTotalHeatFlux);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputeHeatFluxWithNeighbor(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Compute simulated or adjusted interaction properties
    ComputeInteractionProps(r_process_info);

    // Heat generation
    // ASSUMPTION: Heat is generated even when neighbor is adiabatic
    if (r_process_info[HEAT_GENERATION_OPTION] && mHasMotion)
      mGenerationHeatFlux += GetGenerationModel().ComputeHeatGeneration(r_process_info, this);

    // Check if neighbor is adiabatic
    if (CheckAdiabaticNeighbor())
      return;

    // Heat transfer mechanisms
    if (r_process_info[DIRECT_CONDUCTION_OPTION])
      mConductionDirectHeatFlux += GetDirectConductionModel().ComputeHeatFlux(r_process_info, this);

    if (r_process_info[INDIRECT_CONDUCTION_OPTION])
      mConductionIndirectHeatFlux += GetIndirectConductionModel().ComputeHeatFlux(r_process_info, this);

    if (r_process_info[RADIATION_OPTION])
      mRadiationHeatFlux += GetRadiationModel().ComputeHeatFlux(r_process_info, this);

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
    if (r_process_info[REAL_CONTACT_OPTION] && mNeighborInContact) {
      GetRealContactModel().AdjustContact(r_process_info, this);
    }
    else {
      mNeighborDistanceAdjusted   = mNeighborDistance;
      mNeighborSeparationAdjusted = mNeighborSeparation;
      mContactRadiusAdjusted      = mContactRadius;
    }

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::StoreBallToBallContactInfo(const ProcessInfo& r_process_info,
                                                          SphericParticle::ParticleDataBuffer& data_buffer,
                                                          double GlobalContactForceTotal[3],
                                                          double LocalContactForceTotal[3],
                                                          double LocalContactForceDamping[3],
                                                          bool   sliding) {
    KRATOS_TRY

    SphericParticle::StoreBallToBallContactInfo(r_process_info, data_buffer, GlobalContactForceTotal, LocalContactForceTotal, LocalContactForceDamping, sliding);

    if (!mStoreContactParam)
      return;

    // Increment number of contact particle neighbors
    SphericParticle* neighbor = data_buffer.mpOtherParticle;
    mNumberOfContactParticleNeighbor++;

    // New contact
    if (!mContactParamsParticle.count(neighbor)) {
      // Add new parameters to map
      ContactParams params;
      mContactParamsParticle[neighbor] = params;

      // Initialize contact parameters
      mContactParamsParticle[neighbor].impact_time         = r_process_info[TIME];
      mContactParamsParticle[neighbor].impact_velocity     = { data_buffer.mLocalRelVel[2], sqrt(data_buffer.mLocalRelVel[0] * data_buffer.mLocalRelVel[0] + data_buffer.mLocalRelVel[1] * data_buffer.mLocalRelVel[1]) }; // local components of relavive velocity (normal and tangential)
      mContactParamsParticle[neighbor].viscodamping_energy = 0.0;
      mContactParamsParticle[neighbor].frictional_energy   = 0.0;
      mContactParamsParticle[neighbor].rollresist_energy   = 0.0;
    }

    if (r_process_info[HEAT_GENERATION_OPTION]) {
      // If thermal problem was solved in previous step, reset dissipated energies accumulated for this interaction
      if ((r_process_info[TIME_STEPS] - 2) % r_process_info[THERMAL_FREQUENCY] == 0) {
        mContactParamsParticle[neighbor].viscodamping_energy = 0.0;
        mContactParamsParticle[neighbor].frictional_energy   = 0.0;
        mContactParamsParticle[neighbor].rollresist_energy   = 0.0;
      }

      // Update energy dissipated in this interaction (since last thermal solution) with the difference between
      // current and previous accumulated dissipations
      // ATTENTION: Energy increment is multiplied by the inverse of the partition factor used during energy calculation (0.5 to each particle)
      if (r_process_info[GENERATION_DAMPING_OPTION]) {
        const double current_viscodamping_energy   = GetInelasticViscodampingEnergy();
        const double increment_viscodamping_energy = 2.0 * (current_viscodamping_energy - mPreviousViscodampingEnergy);
        mPreviousViscodampingEnergy                = current_viscodamping_energy;
        mContactParamsParticle[neighbor].viscodamping_energy += increment_viscodamping_energy;
      }

      if (r_process_info[GENERATION_SLIDING_OPTION]) {
        const double current_frictional_energy   = GetInelasticFrictionalEnergy();
        const double increment_frictional_energy = 2.0 * (current_frictional_energy - mPreviousFrictionalEnergy);
        mPreviousFrictionalEnergy                = current_frictional_energy;
        mContactParamsParticle[neighbor].frictional_energy += increment_frictional_energy;
      }

      if (r_process_info[GENERATION_ROLLING_OPTION] && this->Is(DEMFlags::HAS_ROTATION) && this->Is(DEMFlags::HAS_ROLLING_FRICTION)) {
        const double current_rollingresistance_energy   = GetInelasticRollingResistanceEnergy();
        const double increment_rollingresistance_energy = 2.0 * (current_rollingresistance_energy - mPreviousRollResistEnergy);
        mPreviousRollResistEnergy                       = current_rollingresistance_energy;
        mContactParamsParticle[neighbor].rollresist_energy += increment_rollingresistance_energy;
      }
    }

    // Update time step
    mContactParamsParticle[neighbor].updated_step = r_process_info[TIME_STEPS];

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::StoreBallToRigidFaceContactInfo(const ProcessInfo& r_process_info,
                                                               SphericParticle::ParticleDataBuffer& data_buffer,
                                                               double GlobalContactForceTotal[3],
                                                               double LocalContactForceTotal[3],
                                                               double LocalContactForceDamping[3],
                                                               bool   sliding) {
    KRATOS_TRY

    SphericParticle::StoreBallToRigidFaceContactInfo(r_process_info, data_buffer, GlobalContactForceTotal, LocalContactForceTotal, LocalContactForceDamping, sliding);

    if (!mStoreContactParam)
      return;

    // Get neighbor wall
    DEMWall* neighbor = data_buffer.mpOtherRigidFace;

    // New contact
    if (!mContactParamsWall.count(neighbor)) {
      // Add new parameters to map
      ContactParams params;
      mContactParamsWall[neighbor] = params;

      // Initialize contact parameters
      mContactParamsWall[neighbor].impact_time         = r_process_info[TIME];
      mContactParamsWall[neighbor].impact_velocity     = { data_buffer.mLocalRelVel[2], sqrt(data_buffer.mLocalRelVel[0] * data_buffer.mLocalRelVel[0] + data_buffer.mLocalRelVel[1] * data_buffer.mLocalRelVel[1]) }; // local components of relavive velocity (normal and tangential)
      mContactParamsWall[neighbor].viscodamping_energy = 0.0;
      mContactParamsWall[neighbor].frictional_energy   = 0.0;
      mContactParamsWall[neighbor].rollresist_energy   = 0.0;
    }

    if (r_process_info[HEAT_GENERATION_OPTION]) {
      // If thermal problem was solved in previous step, reset dissipated energies accumulated for this interaction
      if ((r_process_info[TIME_STEPS] - 2) % r_process_info[THERMAL_FREQUENCY] == 0) {
        mContactParamsWall[neighbor].viscodamping_energy = 0.0;
        mContactParamsWall[neighbor].frictional_energy   = 0.0;
        mContactParamsWall[neighbor].rollresist_energy   = 0.0;
      }

      // Update energy dissipated in this interaction (since last thermal solution) with the difference between
      // current and previous accumulated dissipations
      // ATTENTION: Energy increment is multiplied by the inverse of the partition factor used during energy calculation (1.0 to the particle)
      if (r_process_info[GENERATION_DAMPING_OPTION]) {
        const double current_viscodamping_energy   = GetInelasticViscodampingEnergy();
        const double increment_viscodamping_energy = current_viscodamping_energy - mPreviousViscodampingEnergy;
        mPreviousViscodampingEnergy                = current_viscodamping_energy;
        mContactParamsWall[neighbor].viscodamping_energy += increment_viscodamping_energy;
      }

      if (r_process_info[GENERATION_SLIDING_OPTION]) {
        const double current_frictional_energy   = GetInelasticFrictionalEnergy();
        const double increment_frictional_energy = current_frictional_energy - mPreviousFrictionalEnergy;
        mPreviousFrictionalEnergy                = current_frictional_energy;
        mContactParamsWall[neighbor].frictional_energy += increment_frictional_energy;
      }

      if (r_process_info[GENERATION_ROLLING_OPTION] && this->Is(DEMFlags::HAS_ROTATION) && this->Is(DEMFlags::HAS_ROLLING_FRICTION)) {
        const double current_rollingresistance_energy   = GetInelasticRollingResistanceEnergy();
        const double increment_rollingresistance_energy = current_rollingresistance_energy - mPreviousRollResistEnergy;
        mPreviousRollResistEnergy                       = current_rollingresistance_energy;
        mContactParamsWall[neighbor].rollresist_energy += increment_rollingresistance_energy;
      }
    }

    // Update time step
    mContactParamsWall[neighbor].updated_step = r_process_info[TIME_STEPS];

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::Move(const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag) {
    // Time integration of motion
    if (mHasMotion)
      SphericParticle::Move(delta_t, rotation_option, force_reduction_factor, StepFlag);

    // Time integration of temperature
    if (mIsTimeToSolve && !mHasFixedTemperature && !this->Is(DEMThermalFlags::IS_ADIABATIC)) {
      GetThermalIntegrationScheme().UpdateTemperature(GetGeometry()[0], delta_t * mNumStepsEval, GetParticleHeatCapacity()); // TODO: Remove last argument (capacity becomes a node property accessed with GetFastProperties - same as MOMENT_OF_INERTIA)
    }
  }

  //=====================================================================================================================================================================================
  // Finalization methods

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::FinalizeSolutionStep(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    if (mHasMotion) {
      SphericParticle::FinalizeSolutionStep(r_process_info);

      // Remove non-contacting neighbors from maps of contact parameters
      if (mStoreContactParam)
        CleanContactParameters(r_process_info);
    }

    // Update temperature dependent radius
    if (mIsTimeToSolve && mHasVariableRadius) {
      double added_search_distance = r_process_info[SEARCH_RADIUS_INCREMENT];
      UpdateTemperatureDependentRadius(r_process_info);
      ComputeAddedSearchDistance(r_process_info, added_search_distance);
      SetSearchRadius(GetParticleRadius() + added_search_distance);
    }

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::UpdateTemperatureDependentRadius(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Update radius
    const double new_radius = GetParticleRadius() * (1.0 + GetParticleExpansionCoefficient() * (GetParticleTemperature() - mPreviousTemperature));
    SetParticleRadius(new_radius);

    // Update inertia
    SetParticleMomentInertia(CalculateMomentOfInertia());

    // TODO: update density

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Auxiliary computations

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::ComputeAddedSearchDistance(const ProcessInfo& r_process_info, double& added_search_distance) {
    KRATOS_TRY

    if (r_process_info[DIRECT_CONDUCTION_OPTION] && mpDirectConductionModel) {
      const double model_search_distance  = GetDirectConductionModel().GetSearchDistance(r_process_info, this);
      const double current_added_distance = added_search_distance;
      added_search_distance = std::max(current_added_distance, model_search_distance);
    }

    if (r_process_info[INDIRECT_CONDUCTION_OPTION] && mpIndirectConductionModel) {
      const double model_search_distance  = GetIndirectConductionModel().GetSearchDistance(r_process_info, this);
      const double current_added_distance = added_search_distance;
      added_search_distance = std::max(current_added_distance, model_search_distance);
    }

    if (r_process_info[RADIATION_OPTION] && mpRadiationModel) {
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

    double distance = 0.0;

    if (mNeighborType & PARTICLE_NEIGHBOR) {
      array_1d<double, 3> direction;
      noalias(direction) = GetParticleCoordinates() - GetNeighborCoordinates();
      distance = DEM_MODULUS_3(direction);
    }
    else if (mNeighborType & WALL_NEIGHBOR_CONTACT) {
      // Computing the distance again, as it is done in SphericParticle::ComputeBallToRigidFaceContactForce
      array_1d<double, 4>& weight = this->mContactConditionWeights[mNeighborIndex];

      // Dummy variables: not used now
      double dummy1[3][3];
      DEM_SET_COMPONENTS_TO_ZERO_3x3(dummy1);
      array_1d<double, 3> dummy2 = ZeroVector(3);
      array_1d<double, 3> dummy3 = ZeroVector(3);
      int dummy4 = 0;

      mNeighbor_w->ComputeConditionRelativeData(mNeighborIndex, this, dummy1, distance, weight, dummy2, dummy3, dummy4);
    }
    else if (mNeighborType & WALL_NEIGHBOR_NONCONTACT) {
      // Computing the distance again, as it is done in SphericParticle::ComputeBallToRigidFaceContactForce

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
    }

    // ATTENTION:
    // If for any reason the distance remain null (some rare cases in particle-wall contact),
    // set it to the summ of the radii of the elements (1% more, as if there is no overlap)
    // to avoid numerical issues of using a zero distance or separation in some formulas.
    if (distance <= std::numeric_limits<double>::epsilon())
      distance = 1.001 * (GetParticleRadius() + GetNeighborRadius()); // GetNeighborRadius should return 0.0 for walls!
    
    return distance;
    
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

    // ATTENTION: Assumption: Original model was not assumed real Young modulus for col_time_max and Rc_max!
    const double col_time_max = ComputeMaxCollisionTimeReal();
    const double Rc_max       = ComputeMaxContactRadiusReal();

    if (Rc_max > 0.0)
      return GetParticleConductivity() * col_time_max / (GetParticleDensity() * GetParticleHeatCapacity() * Rc_max * Rc_max);
    else
      return 0.0;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeMaxCollisionTime(void) {
    KRATOS_TRY

    const double eff_radius             = ComputeEffectiveRadius();
    const double eff_mass               = ComputeEffectiveMass();
    const double eff_young              = ComputeEffectiveYoung();
    const double impact_normal_velocity = fabs(GetContactParameters().impact_velocity[0]);

    if (impact_normal_velocity != 0.0)
      return 2.87 * pow(eff_mass * eff_mass / (eff_radius * eff_young * eff_young * impact_normal_velocity), 0.2);
    else
      return 0.0;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeMaxCollisionTimeReal(void) {
    KRATOS_TRY

    const double eff_radius             = ComputeEffectiveRadius();
    const double eff_mass               = ComputeEffectiveMass();
    const double eff_young              = ComputeEffectiveYoungReal();
    const double impact_normal_velocity = fabs(GetContactParameters().impact_velocity[0]);

    if (impact_normal_velocity != 0.0)
      return 2.87 * pow(eff_mass * eff_mass / (eff_radius * eff_young * eff_young * impact_normal_velocity), 0.2);
    else
      return 0.0;

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
  double ThermalSphericParticle::ComputeMaxContactRadiusReal(void) {
    KRATOS_TRY

    const double eff_radius             = ComputeEffectiveRadius();
    const double eff_mass               = ComputeEffectiveMass();
    const double eff_young              = ComputeEffectiveYoungReal();
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
        if (r >= mNeighborDistance)
            Rc = sqrt(r * r - mNeighborDistance * mNeighborDistance);
        else
            Rc = 0.0;
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

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::ComputeMeanConductivity(void) {
    KRATOS_TRY

    return (GetParticleConductivity() + GetNeighborConductivity()) / 2.0;

    KRATOS_CATCH("")
  }

  //=====================================================================================================================================================================================
  // Get/Set methods

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
  HeatGenerationMechanism& ThermalSphericParticle::GetGenerationModel(void) {
    return *mpGenerationModel;
  }

  //------------------------------------------------------------------------------------------------------------
  RealContactModel& ThermalSphericParticle::GetRealContactModel(void) {
    return *mpRealContactModel;
  }

  //------------------------------------------------------------------------------------------------------------
  ThermalDEMIntegrationScheme& ThermalSphericParticle::GetThermalIntegrationScheme(void) {
    return *mpThermalIntegrationScheme;
  }

  //------------------------------------------------------------------------------------------------------------
  NumericalIntegrationMethod& ThermalSphericParticle::GetNumericalIntegrationMethod(void) {
    return *mpNumericalIntegrationMethod;
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
  array_1d<double, 3> ThermalSphericParticle::GetParticleAngularVelocity(void) {
    return GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
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
  double ThermalSphericParticle::GetParticleCharacteristicLength(void) {
    return 2.0 * GetRadius(); // ATTENTION: What about 2D?
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
  double ThermalSphericParticle::GetParticleDiffusivity(void) {
    return GetParticleConductivity() / (GetParticleDensity() * GetParticleHeatCapacity());
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
  double ThermalSphericParticle::GetContactRollingFrictionCoefficient(void) {
    if (mNeighborType & PARTICLE_NEIGHBOR) {
      Properties& properties_of_contact = GetProperties().GetSubProperties(mNeighbor_p->GetProperties().Id());
      return properties_of_contact[ROLLING_FRICTION];
    }
    else if (mNeighborType & WALL_NEIGHBOR) {
      Properties& properties_of_contact = GetProperties().GetSubProperties(mNeighbor_w->GetProperties().Id());
      return properties_of_contact[ROLLING_FRICTION_WITH_WALLS];
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
      null_param.updated_step        = 0;
      null_param.impact_time         = 0.0;
      null_param.viscodamping_energy = 0.0;
      null_param.frictional_energy   = 0.0;
      null_param.rollresist_energy   = 0.0;
      null_param.impact_velocity.assign(2, 0.0);
      return null_param;
    }
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
  void ThermalSphericParticle::SetGenerationModel(HeatGenerationMechanism::Pointer& model) {
    mpGenerationModel = model->CloneRaw();
  }

  //------------------------------------------------------------------------------------------------------------
  void ThermalSphericParticle::SetRealContactModel(RealContactModel::Pointer& model) {
    mpRealContactModel = model->CloneRaw();
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

  //=====================================================================================================================================================================================
  // DIMENSION DEPENDENT METHODS (DIFFERENT FOR 2D AND 3D)
  // ATTENTION:
  // METHODS INEHERITED IN CYLINDER PARTICLE (2D) FROM SPEHRIC PARTICLE ARE REIMPLEMENTED HERE
  // THIS IS TO AVOID MAKING THERMAL PARTICLE A TEMPALTE CLASS TO INHERIT FROM CYLINDER PARTICLE

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::CalculateVolume(void) {
    const double r = GetParticleRadius();

    if      (mDimension == 2) return Globals::Pi * r * r;
    else if (mDimension == 3) return SphericParticle::CalculateVolume();
    else return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::CalculateMomentOfInertia(void) {
    const double r = GetParticleRadius();

    if      (mDimension == 2) return 0.5 * GetMass() * r * r;
    else if (mDimension == 3) return SphericParticle::CalculateMomentOfInertia();
    else return 0.0;
  }

  //------------------------------------------------------------------------------------------------------------
  double ThermalSphericParticle::GetParticleSurfaceArea(void) {
    const double r = GetParticleRadius();

    if      (mDimension == 2) return 2.0 * Globals::Pi * r;
    else if (mDimension == 3) return 4.0 * Globals::Pi * r * r;
    else return 0.0;
  }

} // namespace Kratos
