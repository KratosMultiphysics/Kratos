//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_AXISYM_WATER_2D_CONDITION_H_INCLUDED)
#define KRATOS_AXISYM_WATER_2D_CONDITION_H_INCLUDED

#include "custom_conditions/axisym_point_rigid_contact_penalty_2D_condition.hpp"

namespace Kratos
{
   class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) AxisymPointRigidContactPenaltyWater2DCondition
      : public  AxisymPointRigidContactPenalty2DCondition
   {
      public:
         typedef Vector VectorType;

         KRATOS_CLASS_POINTER_DEFINITION( AxisymPointRigidContactPenaltyWater2DCondition );

         /// Default constructor.
         AxisymPointRigidContactPenaltyWater2DCondition(IndexType NewId, GeometryType::Pointer pGeometry);

         AxisymPointRigidContactPenaltyWater2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

         AxisymPointRigidContactPenaltyWater2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall);


         /// Copy constructor
         AxisymPointRigidContactPenaltyWater2DCondition( AxisymPointRigidContactPenaltyWater2DCondition const& rOther);


         /// Destructor.
         virtual ~AxisymPointRigidContactPenaltyWater2DCondition();

         /**
          * creates a new condition pointer
          * @param NewId: the ID of the new condition
          * @param ThisNodes: the nodes of the new condition
          * @param pProperties: the properties assigned to the new condition
          * @return a Pointer to the new condition
          */
         Condition::Pointer Create(IndexType NewId, NodesArrayType const&
               ThisNodes,  PropertiesType::Pointer pProperties) const;



         void GetDofList(DofsVectorType& rConditionDofList,
               ProcessInfo& rCurrentProcessInfo );

         /**
          * Sets on rResult the ID's of the element degrees of freedom
          */
         void EquationIdVector(EquationIdVectorType& rResult,
               ProcessInfo& rCurrentProcessInfo );

         /**
          * Sets on rValues the nodal displacements
          */
         void GetValuesVector(Vector& rValues,
               int Step = 0 );

         /**
          * Sets on rValues the nodal velocities
          */
         void GetFirstDerivativesVector(Vector& rValues,
               int Step = 0 );

         /**
          * Sets on rValues the nodal accelerations
          */
         void GetSecondDerivativesVector(Vector& rValues,
               int Step = 0 );

      protected:
         ///@name Protected static Member Variables
         ///@{
         virtual void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
					  VectorType& rRightHandSideVector,
					  Flags& rCalculationFlags);


         // A protected default constructor necessary for serialization
         AxisymPointRigidContactPenaltyWater2DCondition() {};

         virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
               ConditionVariables& rVariables,
               double& rIntegrationWeight);

         virtual void CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix,
               ConditionVariables& rVariables,
               double& rIntegrationWeight);

         /**
          * Calculation of the Contact Forces Vector
          */
         virtual void CalculateAndAddContactForces(Vector& rRightHandSideVector,
               ConditionVariables& rVariables,
               double& rIntegrationWeight );

      private:
         friend class Serializer;

         virtual void save(Serializer& rSerializer) const
         {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, AxisymPointRigidContactPenalty2DCondition )
         }

         virtual void load(Serializer& rSerializer)
         {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, AxisymPointRigidContactPenalty2DCondition )
         }



   }; // end class


   } // end Namespace Kratos



#endif
