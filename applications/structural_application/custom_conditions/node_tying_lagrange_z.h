//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2008-10-23 12:26:35 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_NODE_TYING_LAGRANGE_Z_H_INCLUDED )
#define  KRATOS_NODE_TYING_LAGRANGE_Z_H_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/node.h"
#include "includes/dof.h"
#include "containers/variables_list_data_value_container.h"


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
  class NodeTyingLagrangeZ
	  : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
        
      /// Counted pointer of NodeTyingLagrangeZ
      KRATOS_CLASS_POINTER_DEFINITION(NodeTyingLagrangeZ);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
	  NodeTyingLagrangeZ(IndexType NewId, GeometryType::Pointer pGeometry);
      NodeTyingLagrangeZ(IndexType NewId, GeometryType::Pointer pGeometry,  
PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~NodeTyingLagrangeZ();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      Condition::Pointer Create(IndexType NewId, NodesArrayType const& 
ThisNodes,  PropertiesType::Pointer pProperties) const;

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& 
rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& 
rCurrentProcessInfo);
      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& 
rCurrentProcessInfo);

	  void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& 
CurrentProcessInfo);

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
		
        
	friend class Serializer;

	// A private default constructor necessary for serialization 
	NodeTyingLagrangeZ(){}; 

	virtual void save(Serializer& rSerializer) const
	{
	rSerializer.save("Name","NodeTyingLagrangeZ");
	KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
	}

	virtual void load(Serializer& rSerializer)
	{
	KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
	}
	
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
      
      /// Assignment operator.
      //NodeTyingLagrangeZ& operator=(const NodeTyingLagrangeZ& rOther);

      /// Copy constructor.
      //NodeTyingLagrangeZ(const NodeTyingLagrangeZ& rOther);

        
      ///@}    
        
    }; // Class NodeTyingLagrangeZ 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
				    NodeTyingLagrangeZ& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
				    const NodeTyingLagrangeZ& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_NODE_TYING_LAGRANGE_Z_H_INCLUDED defined 

 

