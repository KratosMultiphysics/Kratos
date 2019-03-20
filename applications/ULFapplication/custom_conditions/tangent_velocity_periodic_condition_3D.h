#if !defined(KRATOS_TANGENT_VELOCITY_LAGRANGE_MULTIPLIER_CONDITION3D_H_INCLUDED )
#define  KRATOS_TANGENT_VELOCITY_LAGRANGE_MULTIPLIER_CONDITION3D_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"


namespace Kratos
{
 class TangentVelocityPeriodicCondition3D2N
  : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
      
       /// Counted pointer of TangentVelocityPeriodicCondition3D2N
       KRATOS_CLASS_POINTER_DEFINITION(TangentVelocityPeriodicCondition3D2N);
      
      /// Default constructor. 
      TangentVelocityPeriodicCondition3D2N(IndexType NewId, GeometryType::Pointer pGeometry);

      TangentVelocityPeriodicCondition3D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~TangentVelocityPeriodicCondition3D2N();
      

      Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateLocalVelocityContribution(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

      void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo);

   protected:
 
 
   private:

	friend class Serializer;

	// A private default constructor necessary for serialization  
	TangentVelocityPeriodicCondition3D2N() : Condition()
       {
       }
        

  }; // Class TangentVelocityPeriodicCondition3D2N 

} //namespace kratos 
#endif  
