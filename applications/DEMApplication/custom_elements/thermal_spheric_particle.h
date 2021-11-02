//
//   Project Name:                     ThermalDEM $
//   Last Modified by:    $Author: Ferran Arrufat $
//   Date:                $Date:    February 2015 $
//   Revision:            $Revision:      1.0.0.0 $
//

#if !defined(KRATOS_THERMAL_SPHERIC_PARTICLE_H_INCLUDED)
#define KRATOS_THERMAL_SPHERIC_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "spheric_particle.h"
#include "spheric_continuum_particle.h"

namespace Kratos
{
template <class TBaseElement>
class KRATOS_API(DEM_APPLICATION) ThermalSphericParticle : public TBaseElement
{
  public:
  
  // Pointer definition of ThermalSphericParticle
  KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ThermalSphericParticle);
  
  typedef GlobalPointersVector<Element> ParticleWeakVectorType;
  typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
  typedef GlobalPointersVector<Element>::iterator ParticleWeakIteratorType;

  typedef Node<3> NodeType;
  typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
  typedef std::size_t IndexType;
  typedef Geometry<Node<3>> GeometryType;
  typedef Properties PropertiesType;

  using TBaseElement::mNeighbourElements;
  using TBaseElement::GetProperties;
  using TBaseElement::GetGeometry;
  using TBaseElement::GetRadius;
  using TBaseElement::GetMass;
  using TBaseElement::SetRadius;

  // Definitions
  #define PARTICLE_NEIGHBOR  1
  #define WALL_NEIGHBOR      2
  #define STEFAN_BOLTZMANN   5.670374419e-8

  struct IntegrandParams
  {
    double x;
    double p1, p2, p3, p4, p5, p6, p7, p8, p9;
  };

  // Constructor
  ThermalSphericParticle():TBaseElement(){};
  ThermalSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry):TBaseElement(NewId, pGeometry){};
  ThermalSphericParticle(IndexType NewId, NodesArrayType const& ThisNodes):TBaseElement(NewId, ThisNodes){};
  ThermalSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):TBaseElement(NewId, pGeometry, pProperties){};

  Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override {
    return Element::Pointer(new ThermalSphericParticle<TBaseElement>(NewId, GetGeometry().Create(ThisNodes), pProperties));
  };

  // Destructor
  virtual ~ThermalSphericParticle();

  // Initialization methods
  void Initialize(const ProcessInfo& r_process_info) override;
  void InitializeSolutionStep(const ProcessInfo& r_process_info) override;
  void InitializeHeatFluxComputation(const ProcessInfo& r_process_info);

  // Calculate right hand side
  void CalculateRightHandSide(const ProcessInfo& r_process_info, double dt, const array_1d<double, 3>& gravity) override;
  void ComputeHeatFluxes(const ProcessInfo& r_process_info);

  // Finalization methods
  void FinalizeSolutionStep(const ProcessInfo& r_process_info) override;

  // Update methods
  void UpdateTemperature(const ProcessInfo& r_process_info);
  void UpdateTemperatureDependentRadius(const ProcessInfo& r_process_info);
  void UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion(const ProcessInfo& r_process_info, double& thermalDeltDisp, double& thermalRelVel, ThermalSphericParticle<TBaseElement>* element2);
  void RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(const ProcessInfo& r_process_info, double DeltDisp[3], double RelVel[3], double OldLocalCoordSystem[3][3], double LocalCoordSystem[3][3], SphericParticle* neighbor_iterator) override;

  // Heat fluxes computation
  void ComputeHeatFluxWithNeighbor(const ProcessInfo& r_process_info);
  void ComputeInteractionProps(const ProcessInfo& r_process_info);
  void ComputeDirectConductionHeatFlux(const ProcessInfo& r_process_info);
  void ComputeIndirectConductionHeatFlux(const ProcessInfo& r_process_info);
  void ComputeRadiativeHeatFlux(const ProcessInfo& r_process_info);
  void ComputeContinuumRadiativeHeatFlux(const ProcessInfo& r_process_info);
  void ComputeFrictionHeatFlux(const ProcessInfo& r_process_info);
  void ComputeConvectiveHeatFlux(const ProcessInfo& r_process_info);
  void ComputePrescribedHeatFlux(const ProcessInfo& r_process_info);

  // Heat transfer models
  double DirectConductionBatchelorOBrien(const ProcessInfo& r_process_info);
  double DirectConductionThermalPipe(const ProcessInfo& r_process_info);
  double DirectConductionCollisional(const ProcessInfo& r_process_info);
  double IndirectConductionSurroundingLayer(const ProcessInfo& r_process_info);
  double IndirectConductionVoronoiA(const ProcessInfo& r_process_info);
  double IndirectConductionVoronoiB(const ProcessInfo& r_process_info);
  double IndirectConductionVargasMcCarthy(const ProcessInfo& r_process_info);
  double NusseltHanzMarshall(const ProcessInfo& r_process_info);
  double NusseltWhitaker(const ProcessInfo& r_process_info);
  double NusseltGunn(const ProcessInfo& r_process_info);
  double NusseltLiMason(const ProcessInfo& r_process_info);
  double RadiationContinuumZhou(const ProcessInfo& r_process_info);
  double RadiationContinuumKrause(const ProcessInfo& r_process_info);
  double AdjustedContactRadiusZhou(const ProcessInfo& r_process_info);
  double AdjustedContactRadiusLu(const ProcessInfo& r_process_info);

  // Auxiliary computations
  void   ComputeAddedSearchDistance(const ProcessInfo& r_process_info, double& added_search_distance) override;
  double ComputePrandtlNumber(const ProcessInfo& r_process_info);
  double ComputeReynoldNumber(const ProcessInfo& r_process_info);
  double ComputeFluidRelativeVelocity(const ProcessInfo& r_process_info);

  // Numerical integration
  double AdaptiveSimpsonIntegration(const ProcessInfo& r_process_info, double a, double b, IntegrandParams params, double (ThermalSphericParticle::*evalIntegrand)(IntegrandParams));
  double RecursiveSimpsonIntegration(double a, double b, double fa, double fb, double fc, double tol, IntegrandParams params, double (ThermalSphericParticle::*evalIntegrand)(IntegrandParams));
  double EvalIntegrandSurrLayer(IntegrandParams params);
  double EvalIntegrandVoronoiWall(IntegrandParams params);
  double EvalIntegrandVoronoiMono(IntegrandParams params);
  double EvalIntegrandVoronoiMulti(IntegrandParams params);

  // Neighbor interaction computations
  virtual void ComputeContactArea(const double rmin, double indentation, double& calculation_area);
  bool   CheckAdiabaticNeighbor();
  bool   CheckSurfaceDistance(const double radius_factor);
  double ComputeDistanceToNeighbor();
  double ComputeDistanceToNeighborAdjusted();
  double ComputeFourierNumber();
  double ComputeMaxCollisionTime();
  double ComputeMaxContactRadius();
  double ComputeContactRadius();
  double ComputeEffectiveRadius();
  double ComputeEffectiveMass();
  double ComputeEffectiveYoung();
  double ComputeEffectiveYoungReal();
  double ComputeEffectiveConductivity();
  double ComputeAverageConductivity();

  // Get/Set methods
  double GetYoung()   override;
  double GetPoisson() override;
  double GetDensity() override;

  array_1d<double,3> GetParticleCoordinates();
  array_1d<double,3> GetParticleVelocity();
  double             GetParticleTemperature();
  double             GetParticleRadius();
  double             GetParticleSurfaceArea();
  double             GetParticleCharacteristicLength();
  double             GetParticleVolume();
  double             GetParticleYoung();
  double             GetParticlePoisson();
  double             GetParticleDensity();
  double             GetParticleMass();
  double             GetParticleHeatCapacity();
  double             GetParticleConductivity();
  double             GetParticleEmissivity();
  double             GetParticleExpansionCoefficient();
  
  array_1d<double,3> GetWallCoordinates();
  double             GetWallTemperature();
  double             GetWallRadius();
  double             GetWallYoung();
  double             GetWallPoisson();
  double             GetWallDensity();
  double             GetWallMass();
  double             GetWallHeatCapacity();
  double             GetWallConductivity();
  double             GetWallEmissivity();

  array_1d<double,3> GetNeighborCoordinates();
  double             GetNeighborTemperature();
  double             GetNeighborRadius();
  double             GetNeighborYoung();
  double             GetNeighborPoisson();
  double             GetNeighborDensity();
  double             GetNeighborMass();
  double             GetNeighborHeatCapacity();
  double             GetNeighborConductivity();
  double             GetNeighborEmissivity();
  
  double             GetContactDynamicFrictionCoefficient();

  void               SetParticleTemperature               (const double temperature);
  void               SetParticleHeatFlux                  (const double heat_flux);
  void               SetParticlePrescribedHeatFluxSurface (const double heat_flux);
  void               SetParticlePrescribedHeatFluxVolume  (const double heat_flux);
  void               SetParticleRadius                    (const double radius);
  void               SetParticleMomentInertia             (const double moment_inertia);
  void               SetParticleRealYoungRatio            (const double ratio);

  // Turn back information as a string.
  virtual std::string Info() const override {
    std::stringstream buffer;
    buffer << "ThermalSphericParticle";
    return buffer.str();
  }

  // Print object information
  virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "ThermalSphericParticle";}
  virtual void PrintData(std::ostream& rOStream) const override {}

  protected:

  // General
  bool is_time_to_solve;

  // Heat flux components
  double mConductionDirectHeatFlux;
  double mConductionIndirectHeatFlux;
  double mRadiationHeatFlux;
  double mFrictionHeatFlux;
  double mConvectionHeatFlux;
  double mPrescribedHeatFluxSurface;
  double mPrescribedHeatFluxVolume;
  double mPrescribedHeatFlux;
  double mTotalHeatFlux;

  // Neighbor data
  ThermalSphericParticle<TBaseElement>* mNeighbor_p;
  DEMWall*                              mNeighbor_w;
  int                                   mNeighborType;

  // Interaction properties
  bool   mNeighborInContact;         // flag for contact interaction
  double mRealYoungRatio;            // real value of Young modulus
  double mContactRadius;             // simulation contact radius
  double mNeighborDistance;          // simulation neighbor distance
  double mContactRadiusAdjusted;     // adjusted contact radius from real Young modulus
  double mNeighborDistanceAdjusted;  // adjusted neighbor distance from adjusted contact radius

  // Radiation environment-related
  int    mEnvironmentCount;
  double mEnvironmentTemperature;
  double mEnvironmentTempAux;

  // History-dependent properties
  double mPreviousTemperature;

  private:

  friend class Serializer;

  virtual void save(Serializer& rSerializer) const override {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle);
  }

  virtual void load(Serializer& rSerializer) override {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle);
  }
}; // Class ThermalSphericParticle

  /*
  // input stream function
  template<TBaseElement>
  inline std::istream& operator >> (std::istream& rIStream, ThermalSphericParticle<TBaseElement>& rThis) {return rIStream;}

  // output stream function
  template<TBaseElement>
  inline std::ostream& operator << (std::ostream& rOStream, const ThermalSphericParticle<TBaseElement>& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  */

} // namespace Kratos

#endif // KRATOS_THERMAL_SPHERIC_PARTICLE_H_INCLUDED defined
