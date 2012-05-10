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
      SphericParticle() : DiscreteElement() {

       mContinuumGroup = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_CONTINUUM);
       mFailureId = !(mContinuumGroup); // if ContinuumGroup != 0 --> mFailureId = 0; mFailureId is 1 when mContinuumGroup=0;
       
      }
      SphericParticle( IndexType NewId, GeometryType::Pointer pGeometry );
      SphericParticle( IndexType NewId, NodesArrayType const& ThisNodes);
      SphericParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );
      
      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;
         
      /// Destructor.
      virtual ~SphericParticle();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      void Initialize();
      void CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo);
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
      void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
      void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);
      void GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo );
      void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);
      void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);
      void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo);
      







      ///***********************************************////////////// AIXO ES DECLARA AKI O LA INITIALITZACIÓ.

      //l'he de definir aqui aquest???
 
       std::size_t& GetNumberOfNeighbours(){return(GetGeometry()(0)->FastGetSolutionStepValue(NUMBER_OF_NEIGHBOURS));};

       
       //double mInitialDelta;
       //vector<int> mVectorContactFailureId;
       int mContactFailureId;

       //vector< double > mVectorContactInitialDelta; R: cal cridar-ho cada cop per no fer copia!!
       double mContactInitialDelta;

       vector<array_1d<double,3> > mVectorContactForces;
       array_1d<double,3> mContactForces;

       int mContinuumGroup;
       int mFailureId;

       int mSwitch;



       ///***********************************************////////////// AIXO ES DECLARA AKI O LA INITIALITZACIÓ.










      // ParticleWeakVectorType mInitialNeighbours = GetGeometry()[0].GetSolutionStepValue(INITIAL_NEIGHBOUR_ELEMENTS);

      // std::vector<double>   mInitialDelta = GetGeometry()[0].GetSolutionStepValue(PARTICLE_CONTACT_INITIAL_DELTA);
      // std::vector<int>  mContactFailureId = GetGeometry()[0].GetSolutionStepValue(PARTICLE_CONTACT_FAILURE_ID);
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
      
            
      ///@}
      
    protected:

        void ComputeForcesGeneral(const ProcessInfo& rCurrentProcessInfo);
        void SetInitialContacts(int case_opt);
        void norm(double Vector[3]);
        void VectorGlobal2Local(double LocalCoordSystem[3][3], double GlobalVector[3], double LocalVector[3]);
        void VectorLocal2Global(double LocalCoordSystem[3][3], double LocalVector[3], double GlobalVector[3]);
        double DotProduct(double Vector1[3], double Vector2[3]);
        void CrossProduct(double u[3], double v[3], double ReturnVector[3]);
        void ComputeContactLocalCoordSystem(double NormalDirection[3], double LocalCoordSystem[3][3]);
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
        unsigned int mDimension;
        double mDampType;
        double mTimeStep;
        double mInertia;
        double mMomentOfInertia;
        double mRealMass;

        std::vector<double> mForce;

        
       
        
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


