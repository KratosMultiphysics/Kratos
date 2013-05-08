//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Nelson $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_SPHERIC_PARTICLE_H_INCLUDED )
#define  KRATOS_SPHERIC_PARTICLE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "discrete_element.h"



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
  class SphericParticle : public DiscreteElement
    {
    public:

      ///@name Type Definitions
      ///@{

      /// Pointer definition of SphericParticle
      KRATOS_CLASS_POINTER_DEFINITION(SphericParticle);

      typedef WeakPointerVector<Element> ParticleWeakVectorType;  //M: l'he afegit jo.. esta be akesta?
      typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
      typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.

      SphericParticle( IndexType NewId, GeometryType::Pointer pGeometry );
      SphericParticle( IndexType NewId, NodesArrayType const& ThisNodes);
      SphericParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

      /// Destructor.
      virtual ~SphericParticle();


      ///@}
      ///@name Operations
      ///@{
      void Initialize();
      void CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo);
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
      void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
      void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);
      void GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo );
      void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);
      void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo);

      int mContinuumGroup;
      int* mpFailureId;

      ///@}
      ///@name Access
      ///@{


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
          buffer << "SphericParticle" ;
          return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "SphericParticle";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const {}


      ///@}
      ///@name Friends
      ///@{

    protected:

      SphericParticle();

      //void SetInitialContacts(int case_opt, ProcessInfo& rCurrentProcessInfo);
      void ComputeBallToBallContactForce(   array_1d<double, 3>& rContactForce, array_1d<double, 3>& rContactMoment, ProcessInfo& rCurrentProcessInfo);
      void ComputeBallToSurfaceContactForce(array_1d<double, 3>& rContactForce, array_1d<double, 3>& rContactMoment, ProcessInfo& rCurrentProcessInfo);
      void ComputeParticleBlockContactForce(const ProcessInfo& rCurrentProcessInfo);
      void ComputeParticleRotationSpring(const ProcessInfo& rCurrentProcessInfo);
      void ReasignInternalDataPointers();

      //FOR DEM_FEM APP
//
//        void ComputeParticleBlockContactForce_With_Rotation();
//        void ComputeParticleBlockContactForce_Without_Rotation();
//        void FindContactFaceOfBlockForParticle(ParticleWeakIteratorType rObj_2, int & RightFace, double LocalCoordSystem[3][3], double Coeff[4],double &DistPToB);

      //member variables DEM


      int mDampType;
      int mElasticityType;
      int mRotationOption;
      int mRotationDampType;
      int mGlobalVariablesOption;
      int mCriticalTimeOption;
      int mUniformMaterialOption;
      int mInitializedVariablesFlag;
      double mRadius;
      double mRealMass;
      double mMomentOfInertia;
      double mRollingFriction;
      double mYoung;
      double mPoisson;
      double mFrictionAngle;
      double mRestitCoeff;
      double mMagicFactor;
      double mGlobalKn;
      double mGlobalKt;

      //vector<int>&                  m_data_continuum_ini_neighbours_ids;
      //WeakPointerVector<Element >&  m_data_continuum_ini_neighbours;
      WeakPointerVector<Element >&  m_data_neighbours_elements;
      //vector<int>&                  m_data_ini_neighbours_ids;
      VectorArray3Double&           m_data_particle_elastic_contact_forces;


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


      ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{




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

      /*
      /// Assignment operator.
      SphericParticle& operator=(SphericParticle const& rOther)
      {
    return *this;
      }

      /// Copy constructor.
      SphericParticle(SphericParticle const& rOther)
      {
    *this = rOther;
      }
      */

      ///@}

    }; // Class SphericParticle

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
                    SphericParticle& rThis){ return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
                    const SphericParticle& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SPHERIC_PARTICLE_H_INCLUDED  defined


