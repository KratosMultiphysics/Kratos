//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.1 $
//
//
//
#include "custom_utilities/water_pressure_utilities_Jacobian.hpp"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos 
{

   WaterPressureJacobianUtilities::WaterPressureJacobianUtilities()
   {
   }

   // *************** RHS of the Hydromechanical Problem. Add the solid skeleton movement part *****************************
   // **********************************************************************************************************************
   Vector & WaterPressureJacobianUtilities::CalculateMassBalance_AddDisplacementPart( HydroMechanicalVariables & rVariables, VectorType & rRightHandSideVector, const double & rIntegrationWeight)
   {
      KRATOS_TRY

      double ScalingConstant; 
      GetScalingConstant( ScalingConstant, rVariables.GetProperties() );

      // 2. Geometric
      const GeometryType & rGeometry = rVariables.GetGeometry();
      const unsigned int number_of_nodes = rGeometry.PointsNumber();

      double Jacobian_GP = 0;
      const VectorType & rN = rVariables.GetShapeFunctions();
      for (unsigned int i = 0; i < number_of_nodes; i++)
         Jacobian_GP += rN(i) * rGeometry[i].FastGetSolutionStepValue( JACOBIAN );


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            // Delta Jacobian
            const double &  rCurrentJacobian = rGeometry[j].FastGetSolutionStepValue( JACOBIAN);
            const double & rPreviousJacobian = rGeometry[j].FastGetSolutionStepValue( JACOBIAN, 1);
            double DeltaJacobian = rCurrentJacobian - rPreviousJacobian; 
            
            rRightHandSideVector(i) -= rN(i) * rN(j) * (DeltaJacobian / Jacobian_GP) * rIntegrationWeight * ScalingConstant / rVariables.detF0; 

         }

      }

      return rRightHandSideVector; 

      KRATOS_CATCH("")
   }


   // ****** Tangent To Mass conservation, part of the solid skeleton deformation for displ form ********
   // ***************************************************************************************************
   Matrix& WaterPressureJacobianUtilities::ComputeSolidSkeletonDeformationMatrix( HydroMechanicalVariables & rVariables,  MatrixType & rLocalLHS, const double & rIntegrationWeight)
   {
      KRATOS_TRY

      const GeometryType & rGeometry = rVariables.GetGeometry();
      const unsigned int number_of_nodes = rGeometry.PointsNumber();

      // 1. Some constants
      double ScalingConstant;
      GetScalingConstant( ScalingConstant, rVariables.GetProperties() );

      rLocalLHS.resize( number_of_nodes, number_of_nodes);
      noalias( rLocalLHS ) = ZeroMatrix( number_of_nodes, number_of_nodes );

      const VectorType & rN = rVariables.GetShapeFunctions();


      double Jacobian_GP = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
         Jacobian_GP += rN(i) * rGeometry[i].FastGetSolutionStepValue( JACOBIAN );

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < number_of_nodes; j++) {
            rLocalLHS(i,j) += rN(i) * rN(j) / Jacobian_GP; 
         }
      }

      // Some LD Term (?). 
      double DeltaJacobian_GP = 0;

      for ( unsigned int j = 0; j < number_of_nodes; j++) {
            // Delta Jacobian
            const double &  rCurrentJacobian = rGeometry[j].FastGetSolutionStepValue( JACOBIAN);
            const double & rPreviousJacobian = rGeometry[j].FastGetSolutionStepValue( JACOBIAN, 1);
            double DeltaJacobian = rCurrentJacobian - rPreviousJacobian;
            DeltaJacobian_GP += rN(j) * DeltaJacobian; 
      }

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j <number_of_nodes; j++) {
            rLocalLHS(i,j) -= rN(i) * ( DeltaJacobian_GP / pow( Jacobian_GP, 2.0) ) * rN(j);
         }
      }


      rLocalLHS *= ( rIntegrationWeight * ScalingConstant);

      return rLocalLHS;

      KRATOS_CATCH("")
   }

   // ***************** RESHAPE MATRIX: matrix different in Jacobian Elements ******************************************************
   // ******************************************************************************************************************************
   Matrix & WaterPressureJacobianUtilities::AddReshapeSolidSkeletonDeformationMatrix( Matrix & rLeftHandSide, const Matrix & rLocalLHS, const unsigned int dimension, const unsigned int number_of_variables, const unsigned int number_of_nodes)
   {
      KRATOS_TRY

      Matrix Previous = rLeftHandSide ; 
      // K wP J
      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < number_of_nodes; j++) {
            rLeftHandSide( (i+1)*number_of_variables -1, (j+1)*number_of_variables-2 ) += rLocalLHS(i,j);
         }
      }


      return rLeftHandSide; 

      KRATOS_CATCH("")
   }

   double WaterPressureJacobianUtilities::CalculateVolumeChange( const GeometryType & rGeometry, const Vector & rN, const Matrix & rTotalF)
   {
      const unsigned int number_of_nodes = rGeometry.PointsNumber();

      double VolumeChange = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
         VolumeChange += rN(i) * rGeometry[i].FastGetSolutionStepValue( JACOBIAN );

      return VolumeChange; 

   }
}
