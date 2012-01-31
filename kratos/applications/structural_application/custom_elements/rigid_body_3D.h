//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-06-11 11:00:45 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_RIGID_BODY_3D_INCLUDED )
#define  KRATOS_RIGID_BODY_3D_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"


namespace Kratos
{

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
  class RigidBody3D
	  : public Element
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of RigidBody3D
      KRATOS_CLASS_POINTER_DEFINITION(RigidBody3D);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
	  RigidBody3D(IndexType NewId, GeometryType::Pointer pGeometry);

      
      RigidBody3D(IndexType NewId, 
			       GeometryType::Pointer pGeometry,  
	  		       PropertiesType::Pointer pProperties, 
		               GeometryType::Pointer rskin_nodes, 
			       double& mass, 
		               Matrix& Inertia,
				array_1d<double,3>& translational_stiffness,
				array_1d<double,3>& rotational_stiffness);

      /// Destructor.
      virtual ~RigidBody3D();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

	  void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
	  
	  void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo);
	  void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo);

// 	  void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);
// 	  void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);

	  void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);
	  void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
	  void GetValuesVector(Vector& values, int Step = 0);
	  void GetFirstDerivativesVector(Vector& values, int Step = 0);
	  void GetSecondDerivativesVector(Vector& values, int Step = 0);
	  
	  void Initialize();

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
//      virtual String Info() const;
      
      /// Print information about this object.
//      virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;
      
            
      ///@}      
      ///@name Friends
      ///@{

            
      ///@}
      
    protected:
      ///@name Protected static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected member Variables 
      ///@{ 
      void CalculateForces(VectorType& rExtForces, ProcessInfo& rCurrentProcessInfo);
      void UpdateExtShape(ProcessInfo& rCurrentProcessInfo);

        
        
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
     GeometryType::Pointer  m_skin_nodes;
     double mmass;
     Matrix mInertia;
     boost::numeric::ublas::bounded_matrix<double,3,3> mRotatedInertia;		
     boost::numeric::ublas::bounded_matrix<double,3,3> mRot;		
     array_1d<double,3> mtranslational_stiffness;
     array_1d<double,3> mrotational_stiffness;
      
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
	
	private:

	///@}
	///@name Serialization
	///@{	
	friend class Serializer; 

	// A private default constructor necessary for serialization 
	RigidBody3D(){}

	virtual void save(Serializer& rSerializer) const
	{
	  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
	}

	virtual void load(Serializer& rSerializer)
	{
	  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,  Element );
	}
        
      ///@} 
      ///@name Private  Access 
      ///@{ 
        
        
      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      //RigidBody3D& operator=(const RigidBody3D& rOther);

      /// Copy constructor.
      //RigidBody3D(const RigidBody3D& rOther);

        
      ///@}    
        
    }; // Class RigidBody3D 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
				    RigidBody3D& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
				    const RigidBody3D& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_RIGID_BODY_3D_INCLUDED  defined 


