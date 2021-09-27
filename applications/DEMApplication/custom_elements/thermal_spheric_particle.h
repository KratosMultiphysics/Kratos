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
  
  using TBaseElement::GetDensity;
  using TBaseElement::GetGeometry;
  using TBaseElement::GetMass;
  using TBaseElement::GetPoisson;
  using TBaseElement::GetProperties;
  using TBaseElement::GetRadius;
  using TBaseElement::GetValue;
  using TBaseElement::GetYoung;
  using TBaseElement::SetRadius;
  using TBaseElement::SetValue;
  using TBaseElement::mNeighbourElements;

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

  // Calculate right hand side
  void CalculateRightHandSide(const ProcessInfo& r_process_info, double dt, const array_1d<double, 3>& gravity) override;
  void ComputeHeatFluxes(const ProcessInfo& r_process_info);

  // Compute heat fluxes components
  void ComputeBallToBallDirectConductionHeatFlux(const ProcessInfo& r_process_info);
  void ComputeBallToRigidFaceDirectConductionHeatFlux(const ProcessInfo& r_process_info);
  void ComputeBallToBallIndirectConductionHeatFlux(const ProcessInfo& r_process_info);
  void ComputeBallToRigidFaceIndirectConductionHeatFlux(const ProcessInfo& r_process_info);
  void ComputeConvectiveHeatFlux(const ProcessInfo& r_process_info);
  void ComputeRadiativeHeatFlux(const ProcessInfo& r_process_info);

  // Update methods
  void UpdateTemperature(const ProcessInfo& r_process_info);
  void UpdateTemperatureDependentRadius(const ProcessInfo& r_process_info);
  void UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion(const ProcessInfo& r_process_info, double& thermalDeltDisp, double& thermalRelVel, ThermalSphericParticle<TBaseElement>* element2);
  void RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(const ProcessInfo& r_process_info, double DeltDisp[3], double RelVel[3], double OldLocalCoordSystem[3][3], double LocalCoordSystem[3][3], SphericParticle* neighbour_iterator) override;

  // Finalization methods
  void FinalizeSolutionStep(const ProcessInfo& r_process_info) override;

  // Auxiliary computation methods
  virtual void ComputeContactArea(const double rmin, double indentation, double& calculation_area);
  void ComputeAddedSearchDistance(const ProcessInfo& r_process_info, double& added_search_distance) override;
  double IntegralSurrLayer(const ProcessInfo& r_process_info, double a, double b, double r1, double r2, double d);
  double SolveIntegralSurrLayer(const ProcessInfo& r_process_info, double a, double b, double fa, double fb, double fc, double tol, double r1, double r2, double d);
  double EvalIntegrandSurrLayer(const ProcessInfo& r_process_info, double r, double r1, double r2, double d);
  double IntegralVoronoiMono(const ProcessInfo& r_process_info, double a, double b, double rp, double d, double rij, double keff);
  double SolveIntegralVoronoiMono(const ProcessInfo& r_process_info, double a, double b, double fa, double fb, double fc, double tol, double rp, double d, double rij, double keff);
  double EvalIntegrandVoronoiMono(const ProcessInfo& r_process_info, double r, double rp, double d, double rij, double keff);
  double IntegralVoronoiMulti(const ProcessInfo& r_process_info, double a, double b, double r1, double r2, double d, double rij, double rij_, double D1, double D2);
  double SolveIntegralVoronoiMulti(const ProcessInfo& r_process_info, double a, double b, double fa, double fb, double fc, double tol, double r1, double r2, double d, double rij, double rij_, double D1, double D2);
  double EvalIntegrandVoronoiMulti(const ProcessInfo& r_process_info, double r, double r1, double r2, double d, double rij, double rij_, double D1, double D2);

  // Get/Set methods
  const double& GetParticleTemperature();
  void SetParticleTemperature(const double temperature);
  void SetParticlePrescribedHeatFlux(const double heat_flux);

  // Turn back information as a string.
  virtual std::string Info() const override {
    std::stringstream buffer;
    buffer << "ThermalSphericParticle" ;
    return buffer.str();
  }
  
  // Print object information
  virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "ThermalSphericParticle";}
  virtual void PrintData(std::ostream& rOStream) const override {}

  protected:

  double mThermalConductivity;
  double mSpecificHeat;
  double mConductiveHeatFlux;
  double mConvectiveHeatFlux;
  double mRadiativeHeatFlux;
  double mPrescribedHeatFlux;
  double mTotalHeatFlux;
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
