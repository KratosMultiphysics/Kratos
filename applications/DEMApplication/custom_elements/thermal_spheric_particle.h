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
  typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;

  typedef Node <3> NodeType;
  typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
  typedef std::size_t IndexType;
  typedef Geometry<Node < 3 > > GeometryType;
  typedef Properties PropertiesType;
  
  using TBaseElement::GetGeometry;
  using TBaseElement::GetProperties;
  using TBaseElement::mNeighbourElements;
  using TBaseElement::GetRadius;
  using TBaseElement::SetRadius;
  using TBaseElement::GetMass;
  using TBaseElement::GetValue;
  using TBaseElement::SetValue;

  // Default constructor
  ThermalSphericParticle():TBaseElement(){};
  ThermalSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry):TBaseElement(NewId, pGeometry){};
  ThermalSphericParticle(IndexType NewId, NodesArrayType const& ThisNodes):TBaseElement(NewId, ThisNodes){};
  ThermalSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):TBaseElement(NewId, pGeometry, pProperties){};

  Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override {
    return Element::Pointer(new ThermalSphericParticle<TBaseElement>(NewId, GetGeometry().Create(ThisNodes), pProperties));
  };

  // Destructor
  virtual ~ThermalSphericParticle();

  // Get/Set methods
  const double& GetTemperature();
  const double& GetAmbientTemperature();
  void SetTemperature(const double temperature);

  // Initialization methods
  void Initialize(const ProcessInfo& r_process_info) override;
  void InitializeSolutionStep(const ProcessInfo& r_process_info) override;

  // Calculate right hand side
  void CalculateRightHandSide(const ProcessInfo& r_current_process_info, double dt, const array_1d<double, 3>& gravity) override;
  void ComputeHeatFluxes(const ProcessInfo& r_process_info);

  // Compute heat fluxes components
  void ComputeBallToBallDirectConductionHeatFlux(const ProcessInfo& r_process_info);
  void ComputeBallToRigidFaceDirectConductionHeatFlux(const ProcessInfo& r_process_info);
  void ComputeConvectiveHeatFlux(const ProcessInfo& r_process_info);

  // Auxiliary computation methods
  virtual void ComputeContactArea(const double rmin, double indentation, double& calculation_area);

  // Update methods
  void UpdateTemperature(const ProcessInfo& r_process_info);
  void UpdateTemperatureDependentRadius(const ProcessInfo& r_process_info);
  void UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion(const ProcessInfo& r_process_info, double& thermalDeltDisp, double& thermalRelVel, ThermalSphericParticle<TBaseElement>* element2);

  // Finalization methods
  void FinalizeSolutionStep(const ProcessInfo& r_process_info) override;

  // Others
  void RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(const ProcessInfo& r_process_info, double DeltDisp[3], double RelVel[3], double OldLocalCoordSystem[3][3], double LocalCoordSystem[3][3], SphericParticle* neighbour_iterator) override;

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
  double mTotalHeatFlux;
  double mPreviousTemperature;

  private:

  friend class Serializer;

  virtual void save(Serializer& rSerializer) const override {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle );
  }

  virtual void load(Serializer& rSerializer) override {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle );
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
