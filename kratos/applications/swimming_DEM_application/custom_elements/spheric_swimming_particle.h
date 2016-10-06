//
//   Project Name:        Kratos
//   Last Modified by:    $Author: G.Casas, gcasas@cimne.upc.edu $
//   Date:                $Date: 2013-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//
#if !defined(KRATOS_SPHERIC_SWIMMING_PARTICLE_H_INCLUDED)
#define  KRATOS_SPHERIC_SWIMMING_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "../../applications/DEM_application/custom_elements/spheric_particle.h" 
#include "../../applications/DEM_application/custom_elements/nanoparticle.h"

namespace Kratos
{

  template< class TBaseElement >
  class KRATOS_API(SWIMMING_DEM_APPLICATION) SphericSwimmingParticle : public TBaseElement
    {
    public:

      typedef std::size_t IndexType;
      typedef Node <3> NodeType;
      typedef Geometry<NodeType> GeometryType;
      typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
      typedef Properties PropertiesType;
      
      using TBaseElement::GetGeometry;
      using TBaseElement::GetDensity;
      using TBaseElement::mRealMass;
      using TBaseElement::mRadius;
      using TBaseElement::CalculateVolume;
      using TBaseElement::GetMass;
      using TBaseElement::GetForce;

      
      ///@name Type Definitions
      ///@{

      /// Pointer definition of SphericSwimmingParticle
      KRATOS_CLASS_POINTER_DEFINITION(SphericSwimmingParticle);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      SphericSwimmingParticle():TBaseElement(){}
      SphericSwimmingParticle(IndexType NewId, GeometryType::Pointer pGeometry):TBaseElement(NewId, pGeometry){}
      SphericSwimmingParticle(IndexType NewId, NodesArrayType const& ThisNodes):TBaseElement(NewId, ThisNodes){}
      SphericSwimmingParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties):TBaseElement(NewId, pGeometry, pProperties){}

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {          
        return Element::Pointer(new SphericSwimmingParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
      };

      /// Destructor.
      virtual ~SphericSwimmingParticle(){}

      
      void ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force, array_1d<double, 3>& additionally_applied_moment, const ProcessInfo& rCurrentProcessInfo, const array_1d<double,3>& gravity);

      std::vector<Node<3>::Pointer> mNeighbourNodes;
      std::vector<double> mNeighbourNodesDistances;

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
      // variables for Daitche's method
      static std::vector<double> mAjs;
      static std::vector<double> mBns;
      static std::vector<double> mCns;
      static std::vector<double> mDns;
      static std::vector<double> mEns;
      static bool mDaitcheVectorsAreFull;
      // variables for Hinsberg's method
      static double mTimeWindow;
      static std::vector<double> mAs;
      static std::vector<double> mTs;
      static std::vector<double> mAlphas;
      static std::vector<double> mBetas;

      void Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& r_process_info);


    protected:
        
        void ComputeBuoyancy(NodeType& node, array_1d<double, 3>& buoyancy, const array_1d<double,3>& gravity, const ProcessInfo& r_current_process_info);
        void ComputeDragForce(NodeType& node, array_1d<double, 3>& drag_force, const ProcessInfo& r_current_process_info);
        void ComputeVirtualMassForce(NodeType& node, array_1d<double, 3>& virtual_mass_force, const ProcessInfo& r_current_process_info);
        void ComputeBassetForce(NodeType& node, array_1d<double, 3>& basset_force, const ProcessInfo& r_current_process_info);
        void ComputeSaffmanLiftForce(NodeType& node, array_1d<double, 3>& lift_force, const ProcessInfo& r_current_process_info);
        void ComputeMagnusLiftForce(NodeType& node, array_1d<double, 3>& lift_force, const ProcessInfo& r_current_process_info);
        void ComputeHydrodynamicTorque(NodeType& node, array_1d<double, 3>& hydro_torque, const ProcessInfo& r_current_process_info);
        void ComputeBrownianMotionForce(NodeType& node, array_1d<double, 3>& brownian_motion_force, const ProcessInfo& r_current_process_info);
        void ComputeParticleReynoldsNumber(double& r_reynolds);
        void ComputeParticleRotationReynoldsNumber(double r_norm_of_slip_rot, double& r_reynolds);
        void ComputeParticleAccelerationNumber(const array_1d<double, 3>& slip_acc, double& acc_number);
        void MemberDeclarationFirstStep(const ProcessInfo& r_current_process_info);
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
        void UpdateNodalValues(NodeType& node, const array_1d<double, 3>& hydrodynamic_force, const array_1d<double, 3>& hydrodynamic_moment, const array_1d<double, 3>& weight, const array_1d<double, 3>& buoyancy, const array_1d<double, 3>& drag_force, const array_1d<double, 3>& virtual_mass_force, const array_1d<double, 3>& basset_force, const array_1d<double, 3>& saffman_lift_force, const array_1d<double, 3>& magnus_lift_force, const double &force_reduction_coeff, const ProcessInfo& r_current_process_info);
        void ApplyNumericalAveragingWithOldForces(NodeType& node, array_1d<double, 3>& additionally_applied_force, const ProcessInfo& r_current_process_info);
        double ComputeDragCoefficient(const ProcessInfo& r_current_process_info);
        double ComputeStokesDragCoefficient();
        double ComputeWeatherfordDragCoefficient(const ProcessInfo& r_current_process_info);
        void CalculateNewtonianDragCoefficient(int non_newtonian_option, const double reynolds, const double sphericity, double& r_drag_coeff, int drag_modifier_type);
        double CalculateDragCoeffFromSphericity(const double reynolds, double sphericity, int drag_modifier_type);
        double CalculateShahsTerm(double power_law_n,double power_law_k, double power_law_tol, const double& particle_density, double sphericity, int drag_modifier_type);
        double ComputeGanserDragCoefficient(const ProcessInfo& r_current_process_info);
        double ComputeIshiiDragCoefficient(const ProcessInfo& r_current_process_info);
        double ComputeNewtonRegimeDragCoefficient();
        double ComputeIntermediateRegimeDragCoefficient();
        double ComputeHaiderDragCoefficient();
        double ComputeBeetstraDragCoefficient();
        void ComputeGanserParameters(const int isometric_shape, const double dn, double& k_1, double& k_2);
        void ApplyDragPorosityModification(double& drag_coeff);
        double ComputeElSamniLiftCoefficient(const double norm_of_shear_rate, const double vorticity_norm, const ProcessInfo& r_current_process_info);
        double ComputeMeiLiftCoefficient(const double reynolds, const double reynolds_shear);
        double GetDaitcheCoefficient(int order, unsigned int n, unsigned int j, const double last_h_over_h, const int n_steps_per_quad_step);
        void CalculateExplicitFractionalDerivative(NodeType& node, array_1d<double, 3>& fractional_derivative, double& present_coefficient, Vector& historic_integrands, const double last_h_over_h, const int n_steps_per_quad_step);
        void AddHinsbergTailContribution(NodeType& node, array_1d<double, 3>& fractional_derivative_of_slip_vel, const int order, const int n_steps_per_quad_step, const double time, const double delta_time, const double last_h_over_h, vector<double>& historic_integrands);
        void AddHinsbergTailContributionStrict(NodeType& node, array_1d<double, 3>& fractional_derivative_of_slip_vel, const int order, const int n_steps_per_quad_step, const double time, const double delta_time, const double last_h_over_h, vector<double>& historic_integrands);
        double Phi(const double x);
        double Ki(const double alpha, const double beta, const double time);
        void AddFdi(const int order, array_1d<double, 3>& F, const double t_win, const double alpha, const double beta, const double last_h_over_h, const double delta_time, const vector<double>& historic_integrands, const array_1d<double, 3>& oldest_integrand);
        void AddFre(array_1d<double, 3>& old_Fi, const double beta, const double delta_time);
        void Initialize(const ProcessInfo& r_process_info);
        ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{
      
      bool mHasHydroMomentNodalVar;
      bool mHasDragForceNodalVar;
      bool mHasVirtualMassForceNodalVar;
      bool mHasBassetForceNodalVar;
      bool mHasLiftForceNodalVar;
      bool mHasDragCoefficientVar;
      bool mHasOldAdditionalForceVar;
      bool mFirstStep;
      int mCouplingType;
      int mBuoyancyForceType;
      int mDragForceType;
      int mVirtualMassForceType;
      int mBassetForceType;
      int mSaffmanForceType;
      int mMagnusForceType;
      int mFluidModelType;
      int mPorosityCorrectionType;
      int mHydrodynamicTorqueType;
      int mBrownianMotionType;
      int mQuadratureOrder;
      double mFluidDensity;
      double mFluidFraction;
      double mKinematicViscosity;
      double mSphericity;
      double mNormOfSlipVel;
      double mLastTimeStep;
      double mInitialTime;
      double mOldDaitchePresentCoefficient;
      double mLastVirtualMassAddedMass;
      double mLastBassetForceAddedMass;
      array_1d<double, 3> mSlipVel;
      array_1d<double, 3> mOldBassetTerm;

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
  /*inline std::istream& operator >> (std::istream& rIStream, SphericSwimmingParticle& rThis){ return rIStream;}*/

  /// output stream function
  /*inline std::ostream& operator << (std::ostream& rOStream,
                    const SphericSwimmingParticle& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SPHERIC_SWIMMING_PARTICLE_H_INCLUDED  defined


