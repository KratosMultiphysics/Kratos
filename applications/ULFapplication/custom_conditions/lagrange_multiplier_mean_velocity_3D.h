#if !defined(KRATOS_MEAN_VELOCITY_LAGRANGE_MULTIPLIER_CONDITION_3D_H_INCLUDED )
#define  KRATOS_MEAN_VELOCITY_LAGRANGE_MULTIPLIER_CONDITION_3D_H_INCLUDED

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
 class MeanVelocityLagrangeMultiplierCondition3D
  : public Condition
    {
    public:
      ///@name Type Definitions
      ///@{
      
       /// Counted pointer of MeanVelocityLagrangeMultiplierCondition3D
       KRATOS_CLASS_POINTER_DEFINITION(MeanVelocityLagrangeMultiplierCondition3D);
      
      /// Default constructor. 
      MeanVelocityLagrangeMultiplierCondition3D(IndexType NewId, GeometryType::Pointer pGeometry);

      MeanVelocityLagrangeMultiplierCondition3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~MeanVelocityLagrangeMultiplierCondition3D();
      

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
	MeanVelocityLagrangeMultiplierCondition3D() : Condition()
       {
       }
        

  }; // Class MeanVelocityLagrangeMultiplierCondition3D 

} //namespace kratos 
#endif
