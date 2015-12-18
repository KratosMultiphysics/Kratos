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
  class KRATOS_DEMAPPLICATION_EXPORT_DLL SphericSwimmingParticle : public SphericParticle
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
      void ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force, array_1d<double, 3>& additionally_applied_moment, ProcessInfo& rCurrentProcessInfo, const array_1d<double,3>& gravity);

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
        void ComputeBuoyancy(array_1d<double, 3>& buoyancy, const array_1d<double,3>& gravity, ProcessInfo& r_current_process_info);
        void ComputeDragForce(array_1d<double, 3>& drag_force, ProcessInfo& r_current_process_info);
        void ComputeVirtualMassForce(array_1d<double, 3>& virtual_mass_force, ProcessInfo& r_current_process_info);
        void ComputeSaffmanLiftForce(array_1d<double, 3>& lift_force, ProcessInfo& r_current_process_info);
        void ComputeMagnusLiftForce(array_1d<double, 3>& lift_force, ProcessInfo& r_current_process_info);
        void ComputeHydrodynamicTorque(array_1d<double, 3>& hydro_torque, ProcessInfo& r_current_process_info);
        void ComputeParticleReynoldsNumber(double& r_reynolds);
        void ComputeParticleRotationReynoldsNumber(double r_norm_of_slip_rot, double& r_reynolds);
        void ComputeParticleAccelerationNumber(const array_1d<double, 3>& slip_acc, double& acc_number);
        void AdditionalMemberDeclarationFirstStep(const ProcessInfo& r_current_process_info);
        void AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_current_process_info);

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

        void UpdateNodalValues(const array_1d<double, 3>& hydrodynamic_force, const array_1d<double, 3>& hydrodynamic_moment, const array_1d<double, 3>& buoyancy, const array_1d<double, 3>& drag_force, const array_1d<double, 3>& virtual_mass_force, const array_1d<double, 3>& saffman_lift_force, const array_1d<double, 3>& magnus_lift_force);
        double ComputeStokesDragCoefficient();
        double ComputeWeatherfordDragCoefficient(ProcessInfo& r_current_process_info);
        void CalculateNewtonianDragCoefficient(int non_newtonian_option, const double reynolds, const double sphericity, double& r_drag_coeff, int drag_modifier_type);
        double CalculateDragCoeffFromSphericity(const double reynolds, double sphericity, int drag_modifier_type);
        double CalculateShahsTerm(double power_law_n,double power_law_k, double power_law_tol, const double& particle_density, double sphericity, int drag_modifier_type);
        double ComputeGanserDragCoefficient(ProcessInfo& r_current_process_info);
        double ComputeIshiiDragCoefficient(ProcessInfo& r_current_process_info);
        double ComputeNewtonRegimeDragCoefficient();
        double ComputeIntermediateRegimeDragCoefficient();
        double ComputeHaiderDragCoefficient();
        double ComputeBeetstraDragCoefficient();
        void ComputeGanserParameters(const int isometric_shape, const double dn, double& k_1, double& k_2);
        void ApplyDragPorosityModification(double& drag_coeff);
        double ComputeElSamniLiftCoefficient(const double norm_of_shear_rate, const double vorticity_norm, ProcessInfo& r_current_process_info);
        double ComputeMeiLiftCoefficient(const double reynolds, const double reynolds_shear);
        void CustomInitialize();
        ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{

      bool mHasHydroMomentNodalVar;
      bool mHasDragForceNodalVar;
      bool mHasVirtualMassForceNodalVar;
      bool mHasLiftForceNodalVar;
      int mBuoyancyForceType;
      int mDragForceType;
      int mVirtualMassForceType;
      int mSaffmanForceType;
      int mMagnusForceType;
      int mFluidModelType;
      int mPorosityCorrectionType;
      int mHydrodynamicTorqueType;
      double mFluidDensity;
      double mFluidFraction;
      double mKinematicViscosity;
      double mSphericity;
      double mNormOfSlipVel;
      array_1d<double, 3> mSlipVel;

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


