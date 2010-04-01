//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: jmora $
//   Date:                $Date: 2010-02-02  $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_POINTCURRENT2D_CONDITION_H_INCLUDED )
#define  KRATOS_POINTCURRENT2D_CONDITION_H_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
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
  class PointCurrent2D
	  : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of PointCurrent2D
      KRATOS_CLASS_POINTER_DEFINITION(PointCurrent2D);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
	  PointCurrent2D(IndexType NewId, GeometryType::Pointer pGeometry);
      PointCurrent2D(IndexType NewId, GeometryType::Pointer pGeometry,  
PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~PointCurrent2D();
      

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
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      //PointCurrent2D& operator=(const PointCurrent2D& rOther);

      /// Copy constructor.
      //PointCurrent2D(const PointCurrent2D& rOther);

        
      ///@}    
        
    }; // Class PointCurrent2D 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
				    PointCurrent2D& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
				    const PointCurrent2D& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_POINTCURRENT2D_CONDITION_H_INCLUDED  defined 

 

