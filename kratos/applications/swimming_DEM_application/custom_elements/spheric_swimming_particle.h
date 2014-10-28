//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Nelson $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_SPHERIC_SWIMMING_PARTICLE_H_INCLUDED)
#define  KRATOS_SPHERIC_SWIMMING_PARTICLE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "../../applications/DEM_application/custom_elements/spheric_particle.h" 

namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
  */
  class SphericSwimmingParticle : public SphericParticle
    {
    public:


      ///@name Type Definitions
      ///@{

      /// Pointer definition of SphericSwimmingParticle
      KRATOS_CLASS_POINTER_DEFINITION(SphericSwimmingParticle);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.

      SphericSwimmingParticle(IndexType NewId, GeometryType::Pointer pGeometry);
      SphericSwimmingParticle(IndexType NewId, NodesArrayType const& ThisNodes);
      SphericSwimmingParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

      /// Destructor.
      virtual ~SphericSwimmingParticle();


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      ///@}
      ///@name Access
      ///@{
      std::vector<Node<3>::Pointer> mNeighbourNodes;
      std::vector<double>   mNeighbourNodesDistances;

      ///@}
      ///@name Inquiry
      ///@{


      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
        std::stringstream buffer;
        buffer << "SphericSwimmingParticle" ;
        return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "SphericSwimmingParticle";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const {}


      ///@}
      ///@name Friends
      ///@{


      ///@}

    protected:

        SphericSwimmingParticle();
        void ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force, array_1d<double, 3>& additionally_applied_moment, ProcessInfo& rCurrentProcessInfo, const array_1d<double,3>& gravity);
        void ComputeBuoyancy(array_1d<double, 3>& buoyancy, const double& fluid_density, const array_1d<double,3>& gravity, ProcessInfo& rCurrentProcessInfo);
        void ComputeDragForce(array_1d<double, 3>& drag_force, const double& fluid_density, ProcessInfo& rCurrentProcessInfo);
        void ComputeVirtualMassForce(array_1d<double, 3>& virtual_mass_force, const double& fluid_density, ProcessInfo& rCurrentProcessInfo);
        void ComputeLiftForce(array_1d<double, 3>& lift_force, const double& fluid_density, ProcessInfo& rCurrentProcessInfo);
        void ComputeParticleReynoldsNumber(double rNormOfSlipVel, double rViscosity, double& rReynolds);
        void AdditionalMemberDeclarationFirstStep(const ProcessInfo& rCurrentProcessInfo);
        void AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);

      ///@name Protected static Member Variables
      ///@{


      ///@}
      ///@name Protected member Variables
      ///@{


      ///@}
      ///@name Protected Operators
      ///@{


      ///@}
      ///@name Protected Operations
      ///@{


      ///@}
      ///@name Protected  Access
      ///@{


      ///@}
      ///@name Protected Inquiry
      ///@{


      ///@}
      ///@name Protected LifeCycle
      ///@{


      ///@}

    private:

        void UpdateNodalValues(const array_1d<double, 3>& hydrodynamic_force, const array_1d<double, 3>& buoyancy, const array_1d<double, 3>& drag_force, const array_1d<double, 3>& virtual_mass_force, const array_1d<double, 3>& lift_force);
        double ComputeStokesDragCoefficient(const double fluid_density, ProcessInfo& rCurrentProcessInfo);
        double ComputeWeatherfordDragCoefficient(const double& norm_of_slip_vel, const double fluid_density, ProcessInfo& rCurrentProcessInfo);
        void CalculateNewtonianDragCoefficient(int NonNewtonianOption, const double Reynolds, const double Sphericity, double& rDrag_coeff, int DragModifierType);
        double CalculateDragCoeffFromSphericity(const double Reynolds, double Sphericity, int DragModifierType);
        double CalculateShahsTerm(double PowerLawN,double PowerLawK, double PowerLawTol, const double& ParticleDensity, const double& FluidDensity, double Sphericity, int DragModifier);
        double ComputeGanserDragCoefficient(const double& norm_of_slip_vel, const double fluid_density, ProcessInfo& rCurrentProcessInfo);
        double ComputeIshiiDragCoefficient(const double& norm_of_slip_vel, const double fluid_density, ProcessInfo& rCurrentProcessInfo);
        double ComputeNewtonRegimeDragCoefficient(const double& norm_of_slip_vel, const double fluid_density);
        void ComputeGanserParameters(const int isometric_shape, const double sphericity, const double dn, double& k_1, double& k_2);
        double ComputeSaffmanLiftCoefficient(const double& norm_of_slip_vel, const double fluid_density, const double norm_of_shear_rate, const double vorticity_norm, ProcessInfo& rCurrentProcessInfo);
        void CustomInitialize();
        ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{

      bool mHasDragForceNodalVar;
      bool mHasVirtualMassForceNodalVar;
      bool mHasLiftForceNodalVar;
      int mBuoyancyForceType;
      int mDragForceType;
      int mVirtualMassForceType;
      int mLiftForceType;
      int mFluidModelType;

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{


      ///@}
      ///@name Private  Access
      ///@{


      ///@}
      ///@name Private Inquiry
      ///@{


      ///@}
      ///@name Un accessible methods
      ///@{


      ///@}
      ///@name Serialization
      ///@{

      friend class Serializer;

      virtual void save(Serializer& rSerializer) const
      {
          KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DiscreteElement );
      }

      virtual void load(Serializer& rSerializer)
      {
          KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DiscreteElement );
      }



      /// Assignment operator.

      /// Copy constructor.


      ///@}

    }; // Class SphericSwimmingParticle

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
                    SphericSwimmingParticle& rThis){ return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
                    const SphericSwimmingParticle& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SPHERIC_SWIMMING_PARTICLE_H_INCLUDED  defined


