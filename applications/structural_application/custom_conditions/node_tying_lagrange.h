//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2008-10-23 12:26:35 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_NODE_TYING_LAGRANGE_H_INCLUDED )
#define  KRATOS_NODE_TYING_LAGRANGE_H_INCLUDED



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
  class NodeTyingLagrange
	  : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
        
      /// Counted pointer of NodeTyingLagrange
      KRATOS_CLASS_POINTER_DEFINITION(NodeTyingLagrange);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
	  NodeTyingLagrange(IndexType NewId, GeometryType::Pointer pGeometry);
      NodeTyingLagrange(IndexType NewId, GeometryType::Pointer pGeometry,  
PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~NodeTyingLagrange();
      

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
      ///@name Serialization
      ///@{	
      friend class Serializer;

      // A private default constructor necessary for serialization 
      NodeTyingLagrange(){}; 

      virtual void save(Serializer& rSerializer)
      {
      rSerializer.save("Name","NodeTyingLagrange");
      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
      }

      virtual void load(Serializer& rSerializer)
      {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
      }
        
      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      //NodeTyingLagrange& operator=(const NodeTyingLagrange& rOther);

      /// Copy constructor.
      //NodeTyingLagrange(const NodeTyingLagrange& rOther);

        
      ///@}    
        
    }; // Class NodeTyingLagrange 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
				    NodeTyingLagrange& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
				    const NodeTyingLagrange& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_NODE_TYING_LAGRANGE_H_INCLUDED defined 

 

