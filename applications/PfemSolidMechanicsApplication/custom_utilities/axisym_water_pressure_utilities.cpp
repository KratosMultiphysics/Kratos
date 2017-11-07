//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.1 $
//
//
//
#include "custom_utilities/axisym_water_pressure_utilities.hpp"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos 
{

   AxisymWaterPressureUtilities::AxisymWaterPressureUtilities()
   {
   }



   // **** Geometric like term due to the variation of the density effect on the linear momentum balance ******
   // *********************************************************************************************************
   Matrix & AxisymWaterPressureUtilities::ComputeDensityChangeTerm( HydroMechanicalVariables & rVariables, Matrix & rLocalLHS, const double & rIntegrationWeight)
   {
      KRATOS_TRY

      const GeometryType & rGeometry = rVariables.GetGeometry();
      const unsigned int number_of_nodes = rGeometry.PointsNumber();
      const unsigned int dimension = rGeometry.WorkingSpaceDimension();
      Vector VolumeForce = rVariables.GetVolumeForce();
      VolumeForce *= rVariables.detF0; // due to the volumechange

      rLocalLHS.resize( dimension*number_of_nodes, dimension*number_of_nodes);
      noalias( rLocalLHS ) = ZeroMatrix( dimension*number_of_nodes, dimension*number_of_nodes);
      double density_mixture0 = rVariables.GetProperties().GetValue(DENSITY);
      if ( density_mixture0 > 0) {
         VolumeForce /= density_mixture0; 
      }
      else {
         return rLocalLHS; 
      }

      double density_water = rVariables.GetProperties().GetValue(DENSITY_WATER);

      const MatrixType & rDN_DX = rVariables.GetShapeFunctionsDerivatives();
      const VectorType & rN = rVariables.GetShapeFunctions();

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int iDim = 0; iDim < dimension; iDim++) {
            for (unsigned int p = 0; p < number_of_nodes; p++) {
               for (unsigned int pDim = 0 ; pDim < dimension; pDim++) {
                  rLocalLHS( i*dimension + iDim, p*dimension + pDim) -= rN(i) * VolumeForce(iDim) * rDN_DX(p, pDim);
                  if ( pDim ==0)
                     rLocalLHS(i*dimension + iDim, p*dimension + pDim) -= rN(i) * VolumeForce(iDim) * rN(p) * (1.0/ rVariables.CurrentRadius);
               }
            }
         }
      }

      rLocalLHS *= (rIntegrationWeight * density_water); 


      return rLocalLHS; 

      KRATOS_CATCH("")
   }


   // ***************** Tangent To water Contribution to internal forces *********************************
   // ****************************************************************************************************
   Matrix & AxisymWaterPressureUtilities::ComputeWaterPressureKUwP( HydroMechanicalVariables & rVariables, MatrixType & rLocalLHS, const double & rIntegrationWeight)
   {
      const unsigned int number_of_nodes = rVariables.GetGeometry().PointsNumber();
      const unsigned int dimension = rVariables.GetGeometry().WorkingSpaceDimension();

      rLocalLHS = ZeroMatrix( number_of_nodes*dimension, number_of_nodes);

      const MatrixType & rDN_DX = rVariables.GetShapeFunctionsDerivatives();
      const VectorType & rN     = rVariables.GetShapeFunctions();

      for (unsigned int i = 0; i < number_of_nodes; i++) {
         for (unsigned int j = 0; j < number_of_nodes; j++) {
            for (unsigned int dim = 0; dim < dimension; dim++) {
               rLocalLHS(i*dimension + dim, j) += rDN_DX(i, dim) * rN(j) * rIntegrationWeight;
               if ( dim == 0 )
                  rLocalLHS(i*dimension+dim, j) += rN(i) * rN(j) * ( 1.0 / rVariables.CurrentRadius) * rIntegrationWeight;
            }
         }
      }

      return rLocalLHS;
   }

   // ****** Tanget To Mass conservation, part of the solid skeleton deformation for displ form *********
   // ***************************************************************************************************
   Matrix & AxisymWaterPressureUtilities::ComputeSolidSkeletonDeformationMatrix( HydroMechanicalVariables & rVariables, Matrix & rLocalLHS, const double & rIntegrationWeight)
   {
      KRATOS_TRY


      const GeometryType & rGeometry = rVariables.GetGeometry();
      const unsigned int number_of_nodes = rGeometry.PointsNumber();
      const unsigned int dimension = rGeometry.WorkingSpaceDimension();

      // 1. Some constants
      double ScalingConstant;
      GetScalingConstant( ScalingConstant, rVariables.GetProperties() );

      rLocalLHS.resize(number_of_nodes, number_of_nodes*dimension, false);
      noalias( rLocalLHS ) = ZeroMatrix( number_of_nodes, number_of_nodes * dimension);

      // VelocityGradient
      const MatrixType & rDN_DX = rVariables.GetShapeFunctionsDerivatives();
      const VectorType & rN     = rVariables.GetShapeFunctions();

      Matrix Velocity_DX = ZeroMatrix(dimension, dimension);
      for (unsigned int k = 0; k < number_of_nodes; k++)
      {
         const array_1d<double, 3 > &  CurrentDisplacement = rGeometry[k].FastGetSolutionStepValue( DISPLACEMENT );
         const array_1d<double, 3 > & PreviousDisplacement = rGeometry[k].FastGetSolutionStepValue( DISPLACEMENT , 1);
         array_1d<double, 3 > DeltaDisplacement = CurrentDisplacement - PreviousDisplacement; 
         for (unsigned int j = 0; j < dimension; j++)
         {
            for (unsigned int i = 0; i < dimension; i++)
            {
               Velocity_DX(i,j) += DeltaDisplacement[i] * rDN_DX(k,j);
            }
         }
      }

      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         for (unsigned int p = 0; p < number_of_nodes; p++) {
            for (unsigned int pDim = 0; pDim < dimension; pDim++) {
               rLocalLHS( i, p*dimension + pDim) += rN(i)*rDN_DX(p, pDim);

               if ( pDim == 0)
                  rLocalLHS(i, p*dimension + pDim) += rN(i) * rN(p) * ( 1.0 / rVariables.CurrentRadius);
            }
         }
      }

      rLocalLHS *= ( rIntegrationWeight * ScalingConstant);
      return rLocalLHS;
      KRATOS_CATCH("")
   }



   // *************** RHS of the Hydromechanical Problem. Add the solid skeleton movement part *****************************
   // **********************************************************************************************************************
   Vector & AxisymWaterPressureUtilities::CalculateMassBalance_AddDisplacementPart( HydroMechanicalVariables & rVariables, VectorType & rRightHandSideVector, const double & rIntegrationWeight)
   {
      KRATOS_TRY

      double ScalingConstant; 
      GetScalingConstant( ScalingConstant, rVariables.GetProperties());

      // 2. Geometric
      const GeometryType & rGeometry = rVariables.GetGeometry();
      const unsigned int number_of_nodes = rGeometry.PointsNumber();
      unsigned int dimension = rGeometry.WorkingSpaceDimension();

      const MatrixType & rDN_DX = rVariables.GetShapeFunctionsDerivatives();
      const VectorType & rN     = rVariables.GetShapeFunctions();

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            // DeltaDisplacement
            const array_1d<double, 3 > & CurrentDisplacement  = rGeometry[j].FastGetSolutionStepValue( DISPLACEMENT );
            const array_1d<double, 3 > & PreviousDisplacement = rGeometry[j].FastGetSolutionStepValue( DISPLACEMENT , 1 );
            array_1d<double, 3 > DeltaDisplacement      = CurrentDisplacement-PreviousDisplacement;

            for ( unsigned int p = 0; p < dimension; ++p )
            {
               // Solid Skeleton volumetric deformation
               rRightHandSideVector(i) -= rN(i)*DeltaDisplacement[p] * rDN_DX(j,p) * rIntegrationWeight * ScalingConstant / rVariables.detF0;
               if (p == 0)
                  rRightHandSideVector(i) -= rN(i) * DeltaDisplacement[p] * rN(j) * (1.0 / rVariables.CurrentRadius) * rIntegrationWeight * ScalingConstant / rVariables.detF0; 

            }

         }

      }

      return rRightHandSideVector; 


      KRATOS_CATCH("")
   }

   void AxisymWaterPressureUtilities::GetVoigtSize(const unsigned int & dimension, unsigned int & voigtsize, unsigned int & principal_dimension)
   {
      voigtsize = 4; 
      principal_dimension = 3; 
   }

}
// end namespace kratos
