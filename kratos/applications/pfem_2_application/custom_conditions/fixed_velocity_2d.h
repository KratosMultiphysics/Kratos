#if !defined(KRATOS_FIXED_VELOCITY_CONDITION_H_INCLUDED )
#define  KRATOS_FIXED_VELOCITY_CONDITION_H_INCLUDED

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
 class FixedVelocity2D
  : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
      
       /// Counted pointer of FixedVelocity2D
       KRATOS_CLASS_POINTER_DEFINITION(FixedVelocity2D);
      
      /// Default constructor. 
      FixedVelocity2D(IndexType NewId, GeometryType::Pointer pGeometry);

      FixedVelocity2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~FixedVelocity2D();
      

      Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

      void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo);

   protected:
 
 
   private:

	friend class Serializer;

	// A private default constructor necessary for serialization  
	FixedVelocity2D() : Condition()
       {
       }
        

  }; // Class FixedVelocity2D 

} //namespace kratos 
#endif
