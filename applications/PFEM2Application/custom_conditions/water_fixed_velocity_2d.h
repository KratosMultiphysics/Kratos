#if !defined(KRATOS_WATER_FIXED_VELOCITY_CONDITION_H_INCLUDED )
#define  KRATOS_WATER_FIXED_VELOCITY_CONDITION_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"


namespace Kratos
{
 class WaterFixedVelocity2D
  : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
      
       /// Counted pointer of WaterFixedVelocity2D
       KRATOS_CLASS_POINTER_DEFINITION(WaterFixedVelocity2D);
      
      /// Default constructor. 
      WaterFixedVelocity2D(IndexType NewId, GeometryType::Pointer pGeometry);

      WaterFixedVelocity2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~WaterFixedVelocity2D();
      

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
	WaterFixedVelocity2D() : Condition()
       {
       }
        

  }; // Class WaterFixedVelocity2D 

} //namespace kratos 
#endif
