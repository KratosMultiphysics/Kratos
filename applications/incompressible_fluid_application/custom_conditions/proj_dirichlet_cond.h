//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: paolo $
//   Date:                $Date: 2009-05-08 15:38:00 $
//   Revision:            $Revision: 1.0 $
//
//

 
#if !defined(KRATOS_PROJ_DIRICHLET_CONDITION_H_INCLUDED )
#define  KRATOS_PROJ_DIRICHLET_CONDITION_H_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/serializer.h"

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
  class ProjDirichletCond
	  : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of ProjDirichletCond
      KRATOS_CLASS_POINTER_DEFINITION(ProjDirichletCond);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      ProjDirichletCond(IndexType NewId, GeometryType::Pointer pGeometry);
      ProjDirichletCond(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
      ProjDirichletCond(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties, array_1d<double,3> point1, array_1d<double,3> point2, array_1d<double,3> vel1, array_1d<double,3> vel2);

      /// Destructor.
      virtual ~ProjDirichletCond();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

	  void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo);

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
	array_1d<double,3> mPoint1;
	array_1d<double,3> mPoint2;		
        
        array_1d<double,3> mVel1;
	array_1d<double,3> mVel2;	

       
      ///@}
      ///@name Serialization
      ///@{

      friend class Serializer;

      // A private default constructor necessary for serialization
      ProjDirichletCond() : Condition()
      {
      }

      virtual void save(Serializer& rSerializer) const
      {
	 
	  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
      }

      virtual void load(Serializer& rSerializer)
      {
	  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
      }	        

      ///@} 
      ///@name Private Operators
      ///@{ 
        
      
      ///@} 
      ///@name Private Operations
      ///@{ 
        void CalculateN_at_Point(Condition::GeometryType& geom, const double xc, const double yc, array_1d<double,3>& N_at_c);
        
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
      //ProjDirichletCond& operator=(const ProjDirichletCond& rOther);

      /// Copy constructor.
      //ProjDirichletCond(const ProjDirichletCond& rOther);

        
      ///@}    
        
    }; // Class ProjDirichletCond 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
				    ProjDirichletCond& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
				    const ProjDirichletCond& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_PROJ_DIRICHLET_CONDITION_H_INCLUDED  defined 


