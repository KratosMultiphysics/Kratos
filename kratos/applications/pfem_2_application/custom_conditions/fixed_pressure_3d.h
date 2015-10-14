#if !defined(KRATOS_FIXED_PRESSURE_3D_CONDITION_H_INCLUDED )
#define  KRATOS_FIXED_PRESSURE_3D_CONDITION_H_INCLUDED

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
 class FixedPressure3D
  : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
      
       /// Counted pointer of FixedPressure3D
       KRATOS_CLASS_POINTER_DEFINITION(FixedPressure3D);
      
      /// Default constructor. 
      FixedPressure3D(IndexType NewId, GeometryType::Pointer pGeometry);

      FixedPressure3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~FixedPressure3D();
      

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
	FixedPressure3D() : Condition()
       {
       }
        

  }; // Class FixedPressure3D 

} //namespace kratos 
#endif
