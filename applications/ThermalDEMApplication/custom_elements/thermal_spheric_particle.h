//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <limits>

// External includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_elements/spheric_particle.h"

// Project includes
#include "thermal_dem_application_variables.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) ThermalSphericParticle : public SphericParticle
  {
    public:
  
      // Pointer definition
      KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ThermalSphericParticle);
  
      typedef GlobalPointersVector<Element>           ParticleWeakVectorType;
      typedef ParticleWeakVectorType::ptr_iterator    ParticleWeakIteratorType_ptr;
      typedef GlobalPointersVector<Element>::iterator ParticleWeakIteratorType;

      typedef Node                                NodeType;
      typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
      typedef std::size_t                         IndexType;
      typedef Geometry<Node>                      GeometryType;
      typedef Properties                          PropertiesType;

      // Definitions
      #define PARTICLE_NEIGHBOR        4   // binary = 100 (to use in bitwise operations)
      #define WALL_NEIGHBOR            3   // binary = 011 (to use in bitwise operations)
      #define WALL_NEIGHBOR_NONCONTACT 2   // binary = 010 (to use in bitwise operations)
      #define WALL_NEIGHBOR_CONTACT    1   // binary = 001 (to use in bitwise operations)

      struct ContactParams
      {
        int                 updated_step;
        double              impact_time;
        double              viscodamping_energy;
        double              frictional_energy;
        double              rollresist_energy;
        std::vector<double> impact_velocity;
      };

      struct ContactParamsHMS
      {
        array_1d<double,3> normal;
      };

      // Constructor
      ThermalSphericParticle();
      ThermalSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry);
      ThermalSphericParticle(IndexType NewId, NodesArrayType const& ThisNodes);
      ThermalSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

      // Destructor
      virtual ~ThermalSphericParticle();

      // Initialization methods
      void Initialize                    (const ProcessInfo& r_process_info) override;
      void InitializeSolutionStep        (const ProcessInfo& r_process_info) override;
      void InitializeHeatFluxComputation (const ProcessInfo& r_process_info);

      // Computation methods
      void CalculateRightHandSide          (const ProcessInfo& r_process_info, double dt, const array_1d<double, 3>& gravity) override;
      void ComputeHeatFluxes               (const ProcessInfo& r_process_info);
      void ComputeHeatFluxWithNeighbor     (const ProcessInfo& r_process_info);
      void HierarchicalMultiscale          (const ProcessInfo& r_process_info);
      //void StoreContactInfoPP              (SphericParticle::ParticleDataBuffer& data_buffer) override;
      //void StoreContactInfoPW              (SphericParticle::ParticleDataBuffer& data_buffer) override;
      void ComputeInteractionProps         (const ProcessInfo& r_process_info);
      void StoreBallToBallContactInfo      (const ProcessInfo& r_process_info, SphericParticle::ParticleDataBuffer& data_buffer, double GlobalContactForceTotal[3], double LocalContactForceTotal[3], double LocalContactForceDamping[3], bool sliding) override;
      void StoreBallToRigidFaceContactInfo (const ProcessInfo& r_process_info, SphericParticle::ParticleDataBuffer& data_buffer, double GlobalContactForceTotal[3], double LocalContactForceTotal[3], double LocalContactForceDamping[3], bool sliding, double identation) override;
      void Move                            (const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag) override;

      // Finalization methods
      void FinalizeSolutionStep             (const ProcessInfo& r_process_info) override;
      void UpdateDeformationRateRadius      (const ProcessInfo& r_process_info);
      void UpdateTemperatureDependentRadius (const ProcessInfo& r_process_info);

      // Auxiliary computations
      void   ComputeAddedSearchDistance   (const ProcessInfo& r_process_info, double& added_search_distance);
      double ComputePrandtlNumber         (const ProcessInfo& r_process_info);
      double ComputeReynoldNumber         (const ProcessInfo& r_process_info);
      double ComputeFluidRelativeVelocity (const ProcessInfo& r_process_info);
      double GetVoronoiCellFaceRadius     (const ProcessInfo& r_process_info);

      // Neighbor interaction computations
      void   CleanContactParameters              (const ProcessInfo& r_process_info);
      bool   CheckAdiabaticNeighbor              (void);
      bool   CheckSurfaceDistance                (const double radius_factor);
      double ComputeDistanceToNeighbor           (void);
      double ComputeDistanceToNeighborAdjusted   (void);
      double ComputeSeparationToNeighbor         (void);
      double ComputeSeparationToNeighborAdjusted (void);
      double ComputeFourierNumber                (void);
      double ComputeMaxCollisionTime             (void);
      double ComputeMaxCollisionTimeReal         (void);
      double ComputeMaxContactRadius             (void);
      double ComputeMaxContactRadiusReal         (void);
      double ComputeContactRadius                (void);
      double ComputeEffectiveRadius              (void);
      double ComputeEffectiveMass                (void);
      double ComputeEffectiveYoung               (void);
      double ComputeEffectiveYoungReal           (void);
      double ComputeEffectiveConductivity        (void);
      double ComputeAverageConductivity          (void);
      double ComputeMeanConductivity             (void);

      // Get/Set methods
      HeatExchangeMechanism&       GetDirectConductionModel      (void);
      HeatExchangeMechanism&       GetIndirectConductionModel    (void);
      HeatExchangeMechanism&       GetConvectionModel            (void);
      HeatExchangeMechanism&       GetRadiationModel             (void);
      HeatGenerationMechanism&     GetGenerationModel            (void);
      RealContactModel&            GetRealContactModel           (void);
      ThermalDEMIntegrationScheme& GetThermalIntegrationScheme   (void);
      NumericalIntegrationMethod&  GetNumericalIntegrationMethod (void);

      double GetYoung   (void) override;
      double GetPoisson (void) override;
      double GetDensity (void) override;

      array_1d<double,3> GetParticleCoordinates               (void);
      array_1d<double,3> GetParticleVelocity                  (void);
      array_1d<double,3> GetParticleAngularVelocity           (void);
      double             GetParticleTemperature               (void);
      double             GetParticleRadius                    (void);
      double             GetParticleCharacteristicLength      (void);
      double             GetParticleVolume                    (void);
      double             GetParticleYoung                     (void);
      double             GetParticlePoisson                   (void);
      double             GetParticleDensity                   (void);
      double             GetParticleMass                      (void);
      double             GetParticleHeatCapacity              (void);
      double             GetParticleConductivity              (void);
      double             GetParticleDiffusivity               (void);
      double             GetParticleEmissivity                (void);
      double             GetParticleExpansionCoefficient      (void);
  
      array_1d<double,3> GetWallCoordinates                   (void);
      double             GetWallTemperature                   (void);
      double             GetWallRadius                        (void);
      double             GetWallSurfaceArea                   (void);
      double             GetWallYoung                         (void);
      double             GetWallPoisson                       (void);
      double             GetWallDensity                       (void);
      double             GetWallMass                          (void);
      double             GetWallHeatCapacity                  (void);
      double             GetWallConductivity                  (void);
      double             GetWallEmissivity                    (void);

      array_1d<double,3> GetNeighborCoordinates               (void);
      double             GetNeighborTemperature               (void);
      double             GetNeighborRadius                    (void);
      double             GetNeighborSurfaceArea               (void);
      double             GetNeighborYoung                     (void);
      double             GetNeighborPoisson                   (void);
      double             GetNeighborDensity                   (void);
      double             GetNeighborMass                      (void);
      double             GetNeighborHeatCapacity              (void);
      double             GetNeighborConductivity              (void);
      double             GetNeighborEmissivity                (void);
  
      double             GetContactDynamicFrictionCoefficient (void);
      double             GetContactRollingFrictionCoefficient (void);
      ContactParams      GetContactParameters                 (void);

      void               SetDirectConductionModel             (HeatExchangeMechanism::Pointer& model);
      void               SetIndirectConductionModel           (HeatExchangeMechanism::Pointer& model);
      void               SetConvectionModel                   (HeatExchangeMechanism::Pointer& model);
      void               SetRadiationModel                    (HeatExchangeMechanism::Pointer& model);
      void               SetGenerationModel                   (HeatGenerationMechanism::Pointer& model);
      void               SetRealContactModel                  (RealContactModel::Pointer& model);
      void               SetThermalIntegrationScheme          (ThermalDEMIntegrationScheme::Pointer& scheme);
      void               SetNumericalIntegrationMethod        (NumericalIntegrationMethod::Pointer& method);
      void               SetParticleTemperature               (const double temperature);
      void               SetParticleHeatFlux                  (const double heat_flux);
      void               SetParticlePrescribedHeatFluxSurface (const double heat_flux);
      void               SetParticlePrescribedHeatFluxVolume  (const double heat_flux);
      void               SetParticleRadius                    (const double radius);
      void               SetParticleMass                      (const double mass);
      void               SetParticleMomentInertia             (const double moment_inertia);
      void               SetParticleRealYoungRatio            (const double ratio);
      void               SetParticleDeformationRate           (const double rate);
      void               SetParticleDeformationRateStart      (const double start);
      void               SetParticleDeformationRateStop       (const double stop);

      // DIMENSION DEPENDENT METHODS (DIFFERENT FOR 2D AND 3D)
      // ATTENTION:
      // METHODS INEHERITED IN CYLINDER PARTICLE (2D) FROM SPEHRIC PARTICLE ARE REIMPLEMENTED HERE
      // THIS IS TO AVOID MAKING THERMAL PARTICLE A TEMPALTE CLASS TO INHERIT FROM CYLINDER PARTICLE
      double CalculateVolume          (void) override;
      double CalculateMomentOfInertia (void) override;
      double GetParticleSurfaceArea   (void);

      // Pointers to auxiliary objects
      HeatExchangeMechanism*       mpDirectConductionModel;
      HeatExchangeMechanism*       mpIndirectConductionModel;
      HeatExchangeMechanism*       mpConvectionModel;
      HeatExchangeMechanism*       mpRadiationModel;
      HeatGenerationMechanism*     mpGenerationModel;
      RealContactModel*            mpRealContactModel;
      ThermalDEMIntegrationScheme* mpThermalIntegrationScheme;
      NumericalIntegrationMethod*  mpNumericalIntegrationMethod;

      // General properties
      unsigned int mDimension;           // dimension (2D or 3D)
      unsigned int mNumStepsEval;        // number of steps passed since last thermal evaluation
      double       mInitialTemperature;  // temperature from the beginning of the simulation
      double       mPreviousTemperature; // temperature from the beginning of the step
      bool         mIsTimeToSolve;       // flag to solve thermal problem in current step
      bool         mComputeForces;       // flag to solve mechanical behavior by computing forces
      bool         mComputeMotion;       // flag to solve mechanical behavior by computing motion (and forces, as a pre-requisite)
      bool         mHasFixedTemperature; // flag for constant temperature
      bool         mHasVariableRadius;   // flag for temperature-dependent radius
      bool         mStoreContactParam;   // flag to store contact parameters with neighbors when solving the mechanical problem

      // Heat flux components
      double mConductionDirectHeatFlux;
      double mConductionIndirectHeatFlux;
      double mRadiationHeatFlux;
      double mGenerationHeatFlux;
      double mGenerationHeatFlux_damp_particle;
      double mGenerationHeatFlux_damp_wall;
      double mGenerationHeatFlux_slid_particle;
      double mGenerationHeatFlux_slid_wall;
      double mGenerationHeatFlux_roll_particle;
      double mGenerationHeatFlux_roll_wall;
      double mConvectionHeatFlux;
      double mPrescribedHeatFluxSurface;
      double mPrescribedHeatFluxVolume;
      double mPrescribedHeatFlux;
      double mTotalHeatFlux;

      // Energy properties
      double mPreviousViscodampingEnergy; // accumulated energy dissipation from previous interaction: viscodamping 
      double mPreviousFrictionalEnergy;   // accumulated energy dissipation from previous interaction: frictional
      double mPreviousRollResistEnergy;   // accumulated energy dissipation from previous interaction: rolling resistance
      double mGenerationThermalEnergy_damp_particle;
      double mGenerationThermalEnergy_damp_wall;
      double mGenerationThermalEnergy_slid_particle;
      double mGenerationThermalEnergy_slid_wall;
      double mGenerationThermalEnergy_roll_particle;
      double mGenerationThermalEnergy_roll_wall;

      // Heat maps
      std::vector<std::vector<std::vector<double>>> mHeatMapGenerationDampingPP;  // Local heat map matrix for heat generaion by damping between particle-particle
      std::vector<std::vector<std::vector<double>>> mHeatMapGenerationDampingPW;  // Local heat map matrix for heat generaion by damping between particle-wall
      std::vector<std::vector<std::vector<double>>> mHeatMapGenerationSlidingPP;  // Local heat map matrix for heat generaion by sliding between particle-particle
      std::vector<std::vector<std::vector<double>>> mHeatMapGenerationSlidingPW;  // Local heat map matrix for heat generaion by sliding between particle-wall
      std::vector<std::vector<std::vector<double>>> mHeatMapGenerationRollingPP;  // Local heat map matrix for heat generaion by rolling between particle-particle
      std::vector<std::vector<std::vector<double>>> mHeatMapGenerationRollingPW;  // Local heat map matrix for heat generaion by rolling between particle-wall

      // Interaction properties
      bool   mNeighborInContact;          // flag for contact interaction
      double mRealYoungRatio;             // real value of Young modulus
      double mContactRadius;              // simulation contact radius
      double mNeighborDistance;           // simulation neighbor distance
      double mNeighborSeparation;         // simulation neighbor separation (negative value indicates an indentation)
      double mContactRadiusAdjusted;      // adjusted contact radius from real Young modulus
      double mNeighborDistanceAdjusted;   // adjusted neighbor distance from adjusted contact radius
      double mNeighborSeparationAdjusted; // adjusted neighbor separation (negative value indicates an indentation)

      // Radiation environment-related
      unsigned int mRadiativeNeighbors;
      double       mEnvironmentTemperature;
      double       mEnvironmentTempAux;

      // Deformation rate
      double mDeformationRate;
      double mDeformationRateStart;
      double mDeformationRateStop;

      // Neighboring data
      ThermalSphericParticle*                      mNeighbor_p;
      DEMWall*                                     mNeighbor_w;
      int                                          mNeighborType;
      unsigned int                                 mNeighborIndex;
      unsigned int                                 mNumberOfContactParticleNeighbor;
      std::map<SphericParticle*, ContactParams>    mContactParamsParticle;
      std::map<DEMWall*, ContactParams>            mContactParamsWall;
      std::map<SphericParticle*, ContactParamsHMS> mContactParamsParticleHMS;
      std::map<DEMWall*, ContactParamsHMS>         mContactParamsWallHMS;

      // Tesselation data
      unsigned int         mDelaunayPointListIndex;
      std::map<int,double> mNeighborVoronoiRadius;

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "ThermalSphericParticle";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "ThermalSphericParticle";}
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Serializer
      friend class Serializer;

      virtual void save(Serializer& rSerializer) const override {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle);
      }

      virtual void load(Serializer& rSerializer) override {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle);
      }

  }; // Class ThermalSphericParticle

  // input stream function
  inline std::istream& operator >> (std::istream& rIStream, ThermalSphericParticle& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator << (std::ostream& rOStream, const ThermalSphericParticle& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos
