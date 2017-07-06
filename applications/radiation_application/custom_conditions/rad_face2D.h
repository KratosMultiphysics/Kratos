//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti
//

 
#if !defined(KRATOS_RAD_FACE2D_CONDITION_H_INCLUDED )
#define  KRATOS_RAD_FACE2D_CONDITION_H_INCLUDED 



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/convection_diffusion_settings.h"

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
  class RadFace2D
    : public Condition
  {
  public:
    ///@name Type Definitions
    ///@{
    
    /// Counted pointer of RadFace2D
    KRATOS_CLASS_POINTER_DEFINITION(RadFace2D);
    
    ///@}
    ///@name Life Cycle 
    ///@{ 
      
    /// Default constructor.
    RadFace2D(IndexType NewId, GeometryType::Pointer pGeometry);
    RadFace2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~RadFace2D();
      

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
		
        
        
    ///@} 
    ///@name Private Operators
    ///@{ 
    void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
		      ProcessInfo& rCurrentProcessInfo,
		      bool CalculateStiffnessMatrixFlag,
		      bool CalculateResidualVectorFlag);
        ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
  RadFace2D() : Condition()
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
      //ThermalFace2D& operator=(const ThermalFace2D& rOther);

      /// Copy constructor.
      //ThermalFace2D(const ThermalFace2D& rOther);

        
      ///@}    
        
    }; // Class ThermalFace2D 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
				    ThermalFace2D& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
				    const ThermalFace2D& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_THERMAL_FACE2D_CONDITION_H_INCLUDED   defined 


