#if !defined(KRATOS_SECOND_TANGENT_VELOCITY_LAGRANGE_MULTIPLIER_VERTEX_CONDITION3D_H_INCLUDED )
#define  KRATOS_SECOND_TANGENT_VELOCITY_LAGRANGE_MULTIPLIER_VERTEX_CONDITION3D_H_INCLUDED

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
 class SecondTangentVelocityPeriodicVertexCondition3D2N
  : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
      
       /// Counted pointer of SecondTangentVelocityPeriodicVertexCondition3D2N
       KRATOS_CLASS_POINTER_DEFINITION(SecondTangentVelocityPeriodicVertexCondition3D2N);
      
      /// Default constructor. 
      SecondTangentVelocityPeriodicVertexCondition3D2N(IndexType NewId, GeometryType::Pointer pGeometry);

      SecondTangentVelocityPeriodicVertexCondition3D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~SecondTangentVelocityPeriodicVertexCondition3D2N();
      

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
	SecondTangentVelocityPeriodicVertexCondition3D2N() : Condition()
       {
       }
        

  }; // Class SecondTangentVelocityPeriodicVertexCondition3D2N 

} //namespace kratos 
#endif  
