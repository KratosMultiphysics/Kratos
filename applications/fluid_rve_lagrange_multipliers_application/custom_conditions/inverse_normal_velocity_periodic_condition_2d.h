#if !defined(KRATOS_INVERSE_NORMAL_VELOCITY_LAGRANGE_MULTIPLIER_CONDITION_H_INCLUDED )
#define  KRATOS_INVERSE_NORMAL_VELOCITY_LAGRANGE_MULTIPLIER_CONDITION_H_INCLUDED

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
 class InverseNormalVelocityPeriodicCondition2D2N
  : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
      
       /// Counted pointer of InverseNormalVelocityPeriodicCondition2D2N
       KRATOS_CLASS_POINTER_DEFINITION(InverseNormalVelocityPeriodicCondition2D2N);
      
      /// Default constructor. 
      InverseNormalVelocityPeriodicCondition2D2N(IndexType NewId, GeometryType::Pointer pGeometry);

      InverseNormalVelocityPeriodicCondition2D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~InverseNormalVelocityPeriodicCondition2D2N();
      

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
	InverseNormalVelocityPeriodicCondition2D2N() : Condition()
       {
       }
        

  }; // Class InverseNormalVelocityPeriodicCondition2D2N 

} //namespace kratos 
#endif  
