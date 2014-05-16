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
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../custom_constitutive/DEM_discontinuum_constitutive_law.h"



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

      //Cfeng,RigidFace
      typedef WeakPointerVector<Condition> ConditionWeakVectorType; 
      typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;
	  

      typedef WeakPointerVector<Element> ParticleWeakVectorType; 
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
      void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
      void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo);
      void GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo );
      void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo);

      double GetRadius();
      void SetRadius(double radius);
      double GetSqrtOfRealMass();
      void SetSqrtOfRealMass(double sqrt_of_real_mass);

      double GetYoung();
      void SetYoungFromProperties(double* young);      
      double GetRollingFriction();
      void SetRollingFrictionFromProperties(double* rolling_friction);
      double GetPoisson();
      void SetPoissonFromProperties(double* poisson);
      //double GetFrictionAngle();
      //void SetFrictionAngleFromProperties(double& friction_angle);
      double GetTgOfFrictionAngle();
      void SetTgOfFrictionAngleFromProperties(double* tg_of_friction_angle);
      double GetLnOfRestitCoeff();
      void SetLnOfRestitCoeffFromProperties(double* ln_of_restit_coeff);
      double GetDensity();
      void SetDensityFromProperties(double* density);   
      
      PropertiesProxy* GetFastProperties();
      void SetFastProperties(PropertiesProxy* pProps);

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


      
      std::vector<SphericParticle*> mNeighbourElements;
      std::vector<SphericParticle*> mTempNeighbourElements;
      ///@}
      ///@name Friends
      ///@{

    protected:

      SphericParticle();

      //void SetInitialContacts(int case_opt, ProcessInfo& rCurrentProcessInfo);
	  
      ///Cfeng,RigidFace
      virtual void ComputeNewRigidFaceNeighboursHistoricalData();
      virtual void ComputeBallToRigidFaceContactForce(   array_1d<double, 3>& rContactForce, array_1d<double, 3>& rContactMoment, array_1d<double, 3>& rElasticForce, array_1d<double, 3>& InitialRotaMoment, ProcessInfo& rCurrentProcessInfo);
      virtual void ComputeRigidFaceToMeVelocity(ConditionWeakIteratorType rObj_2, std::size_t ino, double LocalCoordSystem[3][3],double & DistPToB, array_1d<double, 3 > &other_to_me_vel, int & ContactType);
      //////
      
      virtual void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);
      virtual void MemberDeclarationFirstStep(const ProcessInfo& rCurrentProcessInfo);
      virtual void ComputeNewNeighboursHistoricalData();
      void CalculateMaxIndentation(double& rCurrentMaxIndentation, const double& rTolerance);
      void CalculateKineticEnergy(double& rKineticEnergy);
      void CalculateElasticEnergyOfContacts(double& rElasticEnergy);
      void CalculateMomentum(array_1d<double, 3>& rMomentum);
      void CalculateLocalAngularMomentum(array_1d<double, 3>& rAngularMomentum); 
      virtual void ComputeBallToBallContactForce(   array_1d<double, 3>& rContactForce, array_1d<double, 3>& rContactMoment, array_1d<double, 3>& rElasticForce, array_1d<double, 3>& InitialRotaMoment, ProcessInfo& rCurrentProcessInfo); 
      void ComputeBallToSurfaceContactForce(array_1d<double, 3>& rContactForce, array_1d<double, 3>& rContactMoment, array_1d<double, 3>& InitialRotaMoment, int surface_num, ProcessInfo& rCurrentProcessInfo);
      void ComputeBallToCylinderContactForce(array_1d<double, 3>& rContactForce, array_1d<double, 3>& rContactMoment, array_1d<double, 3>& InitialRotaMoment, int cylinder_num, ProcessInfo& rCurrentProcessInfo);     
      //virtual void ComputeParticleBlockContactForce(const ProcessInfo& rCurrentProcessInfo);
      //virtual void ComputeParticleRotationSpring(   const ProcessInfo& rCurrentProcessInfo);

      virtual void CalculateEquivalentConstitutiveParameters(array_1d<double, 3>& other_to_me_vect,
                                                             const double& other_radius,
                                                             const double& radius_sum,
                                                             double& kn,
                                                             double& kt,
                                                             double& equiv_visco_damp_coeff_normal,
                                                             double& equiv_visco_damp_coeff_tangential,
                                                             double& equiv_tg_of_fri_ang,
                                                             //ParticleWeakIteratorType neighbour_iterator);
                                                             SphericParticle* neighbour_iterator);

      virtual void EvaluateDeltaDisplacement(double DeltDisp[3],
                                double RelVel[3],
                                //double NormalDir[3],
                                //double OldNormalDir[3],
                                double LocalCoordSystem[3][3],
                                double OldLocalCoordSystem[3][3],
                                array_1d<double, 3> &other_to_me_vect,
                                const array_1d<double, 3> &vel,
                                const array_1d<double, 3> &delta_displ,
                                //ParticleWeakIteratorType neighbour_iterator,
                                SphericParticle* neighbour_iterator,
                                double& distance);
      
      virtual void NormalForceCalculation(double LocalElasticContactForce[3], double kn, double indentation);

      virtual void TangentialForceCalculation(double LocalElasticContactForce[3], double LocalDeltDisp[3], const double& kt, const double& equiv_tg_of_fri_ang, bool& sliding);

      virtual void DisplacementDueToRotation(double DeltDesp[3],
                                //double OldNormalDir[3],
                                double OldLocalCoordSystem[3][3],
                                const double &other_radius,
                                const double &dt,
                                const array_1d<double, 3> &angl_vel,
                                //ParticleWeakIteratorType neighbour_iterator);
                                SphericParticle* neighbour_iterator);

      virtual void ComputeMoments(double LocalElasticContactForce[3],
                                  double GlobalElasticContactForces[3],
                                  //double InitialRotaMoment[3],
                                  array_1d<double, 3>& rInitialRotaMoment,
                                  double LocalCoordSystem[3][3],
                                  const double &other_radius,
                                  array_1d<double, 3>& rContactMoment,
                                  //ParticleWeakIteratorType neighbour_iterator);
                                  SphericParticle* neighbour_iterator);

      virtual void CustomInitialize();
      
      virtual double GetInitialDelta(int index);
      
      virtual void ComputeAdditionalForces(array_1d<double, 3>& contact_force, array_1d<double, 3>& contact_moment,
                                           array_1d<double, 3>& externally_applied_force, array_1d<double, 3>& externally_applied_moment, ProcessInfo& rCurrentProcessInfo);

      virtual void AddUpForcesAndProject(double LocalCoordSystem[3][3],
                                double LocalContactForce[3],
                                double LocalElasticContactForce[3],
                                double GlobalContactForce[3],
                                double GlobalElasticContactForce[3],
                                double ViscoDampingLocalContactForce[3],
                                //double ViscoDampingGlobalContactForce[3],
                                array_1d<double, 3> &rContactForce,
                                array_1d<double, 3> &rElasticForce,
                                const double &i_neighbour_count);
      
      virtual void AddUpFEMForcesAndProject(double LocalCoordSystem[3][3],
                                double LocalContactForce[3],
                                double LocalElasticContactForce[3],
                                double GlobalContactForce[3],
                                double GlobalElasticContactForce[3],
                                double ViscoDampingLocalContactForce[3],
                                //double ViscoDampingGlobalContactForce[3],
                                array_1d<double, 3> &rContactForce,
                                array_1d<double, 3> &rElasticForce,
                                const double &iRigidFaceNeighbour);
      
   
      virtual void CalculateViscoDamping(double LocalRelVel[3],
                                                  double ViscoDampingLocalContactForce[3],
                                                  double indentation,
                                                  double equiv_visco_damp_coeff_normal,
                                                  double equiv_visco_damp_coeff_tangential,
                                                  bool sliding);

      virtual void AdditionalMemberDeclarationFirstStep(const ProcessInfo& rCurrentProcessInfo);
      virtual void AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);
      
      DEMDiscontinuumConstitutiveLaw::Pointer mDiscontinuumConstitutiveLaw;
      
      const int mParticleId; // (NOT YET ACTIVE!!) Identifies the particle biunivocally if it has been properly created (i.e., a non-repeated NewId is passed to the constructor)
      //int mInitializedVariablesFlag;
      int mDimension;
      int mDampType;
      int mElasticityType;
      double mRadius;
      double mSqrtOfRealMass;

      /*double mRollingFriction;
      double mYoung;
      double mPoisson;
      double mTgOfFrictionAngle;
      double mLnOfRestitCoeff;            */

      PropertiesProxy *mFastProperties;
      
      //vector<int> mOldNeighbourIds;//candidate to std::vector
      //vector< array_1d<double, 3> > mOldNeighbourContactForces;//candidate to std::vector
      
      std::vector<int> mOldNeighbourIds;
      std::vector< array_1d<double, 3> > mOldNeighbourContactForces;      
      std::vector<int> mTempNeighboursIds;
      std::vector<array_1d<double, 3> > mTempNeighboursContactForces;
      
      //pointers:
      //const int *mpTimeStep; 
      //const int *mpActivateSearch;
	  
      //Cfeng,RigidFace
      std::vector<int> mFemOldNeighbourIds;
      std::vector< array_1d<double, 3> >  mFemOldNeighbourContactForces;
      std::vector<int> mFemTempNeighboursIds;
      std::vector<array_1d<double, 3> > mFemTempNeighboursContactForces;

      std::vector<double>               mFemTempNeighboursDelta;
      std::vector<int>                  mFemTempNeighboursMapping;

      //std::vector<double>              mFemNeighbourDelta;
      
      //ParticleWeakVectorType mrNeighbours;

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
          rSerializer.save("mRadius",mRadius);
          rSerializer.save("mSqrtOfRealMass",mSqrtOfRealMass);

          /*rSerializer.save("mRollingFriction",mRollingFriction);
          rSerializer.save("mYoung",mYoung);
          rSerializer.save("mPoisson",mPoisson);
          rSerializer.save("mTgOfFrictionAngle",mTgOfFrictionAngle);
          rSerializer.save("mLnOfRestitCoeff",mLnOfRestitCoeff);  */
          //rSerializer.save("mFastProperties",mFastProperties);
      }

      virtual void load(Serializer& rSerializer)
      {
          KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DiscreteElement );
          rSerializer.load("mRadius",mRadius);
          rSerializer.load("mSqrtOfRealMass",mSqrtOfRealMass);

          /*rSerializer.load("mRollingFriction",mRollingFriction);
          rSerializer.load("mYoung",mYoung);
          rSerializer.load("mPoisson",mPoisson);
          rSerializer.load("mTgOfFrictionAngle",mTgOfFrictionAngle);
          rSerializer.load("mLnOfRestitCoeff",mLnOfRestitCoeff);  */
          //rSerializer.load("mFastProperties",mFastProperties);
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


