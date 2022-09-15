//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(KRATOS_THERMAL_SPHERIC_PARTICLE_H_INCLUDED)
#define KRATOS_THERMAL_SPHERIC_PARTICLE_H_INCLUDED

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

      typedef Node<3>                             NodeType;
      typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
      typedef std::size_t                         IndexType;
      typedef Geometry<Node<3>>                   GeometryType;
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
        double              rolling_resistance;
        std::vector<double> impact_velocity;
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
      void ComputeInteractionProps         (const ProcessInfo& r_process_info);
      void StoreBallToBallContactInfo      (const ProcessInfo& r_process_info, SphericParticle::ParticleDataBuffer& data_buffer, double GlobalContactForceTotal[3], double LocalContactForceTotal[3], double LocalContactForceDamping[3], bool sliding) override;
      void StoreBallToRigidFaceContactInfo (const ProcessInfo& r_process_info, SphericParticle::ParticleDataBuffer& data_buffer, double GlobalContactForceTotal[3], double LocalContactForceTotal[3], double LocalContactForceDamping[3], bool sliding) override;
      void Move                            (const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag) override;

      // Finalization methods
      void FinalizeForceComputation                                         (ParticleDataBuffer& data_buffer) override;
      void FinalizeSolutionStep                                             (const ProcessInfo& r_process_info) override;
      void UpdateTemperatureDependentRadius                                 (const ProcessInfo& r_process_info);
      void UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion (const ProcessInfo& r_process_info, double& thermalDeltDisp, double& thermalRelVel, SphericParticle* element2);
      void RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons   (const ProcessInfo& r_process_info, double DeltDisp[3], double RelVel[3], double OldLocalCoordSystem[3][3], double LocalCoordSystem[3][3], SphericParticle* neighbor) override;

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
      double ComputeMaxContactRadius             (void);
      double ComputeContactRadius                (void);
      double ComputeEffectiveRadius              (void);
      double ComputeEffectiveMass                (void);
      double ComputeEffectiveYoung               (void);
      double ComputeEffectiveYoungReal           (void);
      double ComputeEffectiveConductivity        (void);
      double ComputeAverageConductivity          (void);
      double ComputeMeanConductivity             (void);

      // Get/Set methods
      ThermalDEMIntegrationScheme& GetThermalIntegrationScheme   (void);
      NumericalIntegrationMethod&  GetNumericalIntegrationMethod (void);
      HeatExchangeMechanism&       GetDirectConductionModel      (void);
      HeatExchangeMechanism&       GetIndirectConductionModel    (void);
      HeatExchangeMechanism&       GetConvectionModel            (void);
      HeatExchangeMechanism&       GetRadiationModel             (void);
      HeatGenerationMechanism&     GetGenerationModel            (void);
      RealContactModel&            GetRealContactModel           (void);

      double GetYoung   (void) override;
      double GetPoisson (void) override;
      double GetDensity (void) override;

      array_1d<double,3> GetParticleCoordinates               (void);
      array_1d<double,3> GetParticleVelocity                  (void);
      array_1d<double,3> GetParticleAngularVelocity           (void);
      double             GetParticleTemperature               (void);
      double             GetParticleRadius                    (void);
      double             GetParticleSurfaceArea               (void);
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

      void               SetThermalIntegrationScheme          (ThermalDEMIntegrationScheme::Pointer& scheme);
      void               SetNumericalIntegrationMethod        (NumericalIntegrationMethod::Pointer& method);
      void               SetDirectConductionModel             (HeatExchangeMechanism::Pointer& model);
      void               SetIndirectConductionModel           (HeatExchangeMechanism::Pointer& model);
      void               SetConvectionModel                   (HeatExchangeMechanism::Pointer& model);
      void               SetRadiationModel                    (HeatExchangeMechanism::Pointer& model);
      void               SetGenerationModel                   (HeatGenerationMechanism::Pointer& model);
      void               SetRealContactModel                  (RealContactModel::Pointer& model);
      void               SetParticleTemperature               (const double temperature);
      void               SetParticleHeatFlux                  (const double heat_flux);
      void               SetParticlePrescribedHeatFluxSurface (const double heat_flux);
      void               SetParticlePrescribedHeatFluxVolume  (const double heat_flux);
      void               SetParticleRadius                    (const double radius);
      void               SetParticleMass                      (const double mass);
      void               SetParticleMomentInertia             (const double moment_inertia);
      void               SetParticleRealYoungRatio            (const double ratio);

      // Pointers to auxiliary objects
      ThermalDEMIntegrationScheme* mpThermalIntegrationScheme;
      NumericalIntegrationMethod*  mpNumericalIntegrationMethod;
      HeatExchangeMechanism*       mpDirectConductionModel;
      HeatExchangeMechanism*       mpIndirectConductionModel;
      HeatExchangeMechanism*       mpConvectionModel;
      HeatExchangeMechanism*       mpRadiationModel;
      HeatGenerationMechanism*     mpGenerationModel;
      RealContactModel*            mpRealContactModel;

      // General properties
      unsigned int mNumStepsEval;        // number of steps passed since last thermal evaluation
      double       mPreviousTemperature; // temperature from the beginning of the step
      bool         mIsTimeToSolve;       // flag to solve thermal problem in current step
      bool         mHasMotion;           // flag to solve mechanical behavior (forces and displacements)
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
      double mThermalViscodampingEnergy;  // accumulated thermal energy generated due to viscodamping dissipation
      double mThermalFrictionalEnergy;    // accumulated thermal energy generated due to frictional dissipation
      double mThermalRollResistEnergy;    // accumulated thermal energy generated due to rolling resistance dissipation
      double mPreviousViscodampingEnergy; // accumulated viscodamping energy dissipation from previous interaction
      double mPreviousFrictionalEnergy;   // accumulated frictional energy dissipation from previous interaction
      double mPreviousRollResistEnergy;   // accumulated rolling resistance energy dissipation from previous step
      double mPreviousRollResistCoeff;    // total rolling resistance coefficient from previous step

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

      // Neighboring data
      ThermalSphericParticle*                   mNeighbor_p;
      DEMWall*                                  mNeighbor_w;
      int                                       mNeighborType;
      unsigned int                              mNeighborIndex;
      unsigned int                              mNumberOfContactParticleNeighbor;
      std::map<SphericParticle*, ContactParams> mContactParamsParticle;
      std::map<DEMWall*, ContactParams>         mContactParamsWall;

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

#endif // KRATOS_THERMAL_SPHERIC_PARTICLE_H_INCLUDED defined
