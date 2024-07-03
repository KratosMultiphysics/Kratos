//
//   Project Name:        Kratos
//   Last Modified by:    $Author: G.Casas, gcasas@cimne.upc.edu $
//   Date:                $Date: 2019-1-02 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//
#if !defined(KRATOS_SWIMMING_PARTICLE_H_INCLUDED)
#define  KRATOS_SWIMMING_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_elements/spheric_particle.h"
#include "custom_elements/nanoparticle.h"
#include "custom_constitutive/hydrodynamic_interaction_law.h"
#include "custom_constitutive/power_law_hydrodynamic_interaction_law.h"

namespace Kratos
{

  template< class TBaseElement >
  class KRATOS_API(SWIMMING_DEM_APPLICATION) SwimmingParticle : public TBaseElement
    {
    public:

      typedef std::size_t IndexType;
      typedef Node NodeType;
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

      /// Pointer definition of SwimmingParticle
      KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SwimmingParticle);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      SwimmingParticle():TBaseElement(){}
      SwimmingParticle(IndexType NewId, GeometryType::Pointer pGeometry):TBaseElement(NewId, pGeometry){}
      SwimmingParticle(IndexType NewId, NodesArrayType const& ThisNodes):TBaseElement(NewId, ThisNodes){}
      SwimmingParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
        :TBaseElement(NewId, pGeometry, pProperties){}

      Element::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override
      {
        return Element::Pointer(new SwimmingParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
      };

      /// Destructor.
      virtual ~SwimmingParticle(){}

      SwimmingParticle<TBaseElement>& operator=(const SwimmingParticle<TBaseElement>& rOther);
      virtual void CreateHydrodynamicInteractionLaws(const ProcessInfo& r_process_info);

      void ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force,
                                   array_1d<double, 3>& additionally_applied_moment,
                                   const ProcessInfo& rCurrentProcessInfo,
                                   const array_1d<double,3>& gravity) override;

      void FinalizeSolutionStep(const ProcessInfo& r_current_process_info) override;

      void UpdateGentleCouplingInitiationCoefficients(NodeType& node,
                                                      const ProcessInfo& r_current_process_info);

      std::vector<Node::Pointer> mNeighbourNodes;
      std::vector<double> mNeighbourNodesDistances;
      /// NEW IMPLEMENTATION
      //std::vector<Element::Pointer> mNeighbourElements;
      std::vector<Element*> mNeighbourElements;
      ///

      /// Turn back information as a string.
      virtual std::string Info() const override
      {
        std::stringstream buffer;
        buffer << "Swimming version of " << TBaseElement::Info();
        return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const  override
      {
          rOStream << "Swimming version of " << TBaseElement::Info();
      }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const  override{}

      void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                     array_1d<double, 3 > & Output,
                     const ProcessInfo& r_current_process_info) override;

    protected:

        void MemberDeclarationFirstStep(const ProcessInfo& r_current_process_info) override;
        void AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_current_process_info) override;
        array_1d<double,3> ComputeWeight(const array_1d<double,3>& gravity, const ProcessInfo& r_process_info) override;
        void AddCentrifugalForces(array_1d<double,3>& weight, const ProcessInfo& r_process_info);
        void AddCoriolisForces(array_1d<double,3>& weight, const ProcessInfo& r_process_info);
        void AddRelativeAccelerationForces(array_1d<double,3>& weight, const ProcessInfo& r_process_info);
        void AddEulerForces(array_1d<double,3>& weight, const ProcessInfo& r_process_info);
        virtual double GetFluidMass();

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
        double CalculateDragCoeffFromSphericity(const double reynolds, double sphericity, int drag_modifier_type);
        void UpdateNodalValues(NodeType& node,
                               const array_1d<double, 3>& hydrodynamic_force,
                               const array_1d<double, 3>& hydrodynamic_moment,
                               const array_1d<double, 3>& weight,
                               const array_1d<double, 3>& buoyancy,
                               const array_1d<double, 3>& drag_force,
                               const array_1d<double, 3>& inviscid_force,
                               const array_1d<double, 3>& undisturbed_flow_force,
                               const array_1d<double, 3>& history_force,
                               const array_1d<double, 3>& vorticity_induced_lift,
                               const array_1d<double, 3>& rotation_induced_lift,
                               const double &force_reduction_coeff,
                               const ProcessInfo& r_current_process_info);
        void ApplyNumericalAveragingWithOldForces(NodeType& node,
                                                  array_1d<double, 3>& additionally_applied_force,
                                                  const ProcessInfo& r_current_process_info);
        void ApplyDragPorosityModification(double& drag_coeff);
        void Initialize(const ProcessInfo& r_process_info) override;
        ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{

      bool mFirstStep;
      bool mUpdateImpulse = true;
      int mPorosityCorrectionType;
      double mFluidDensity;
      double mKinematicViscosity;
      double mSphericity;
      double mNormOfSlipVel;
      array_1d<double, 3> mPseudoSlipVel;
      array_1d<double, 3> mLastForce;
      HydrodynamicInteractionLaw::Pointer mHydrodynamicInteractionLaw;

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

      virtual void save(Serializer& rSerializer) const override
      {
          KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DiscreteElement );
      }

      virtual void load(Serializer& rSerializer) override
      {
          KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DiscreteElement );
      }



      /// Assignment operator.

      /// Copy constructor.


      ///@}

    }; // Class SwimmingParticle

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

}  // namespace Kratos.

#endif // KRATOS_SWIMMING_PARTICLE_H_INCLUDED  defined


