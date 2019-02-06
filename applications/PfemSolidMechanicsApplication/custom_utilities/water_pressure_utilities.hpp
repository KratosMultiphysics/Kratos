//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.1 $
//
//

#if !defined(KRATOS_WATER_PRESSURE_UTILITIES)
#define KRATOS_WATER_PRESSURE_UTILITIES


// System includes

// External includes
//#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "geometries/geometry.h"
#include "includes/properties.h"

namespace Kratos
{


   class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) WaterPressureUtilities
   {

      public:
         typedef Matrix MatrixType;

         typedef Vector VectorType;

         typedef unsigned int IndexType;

         typedef unsigned int SizeType;

         typedef Properties PropertiesType; 

         typedef Node < 3 > NodeType; 
         typedef Geometry <NodeType>   GeometryType; 


         struct HydroMechanicalVariables {

            private:
               const GeometryType*   mpGeometry;
               const PropertiesType* mpProperties;

               const MatrixType* mpB;
               const MatrixType* mpF0;
               const MatrixType* mpDN_DX;
               const VectorType* mpN;
               const VectorType* mpVolumeForce;

            public:
               double DeltaTime;
               double detF0;
               double CurrentRadius;
               double ConstrainedModulus;

               unsigned int number_of_variables;

               HydroMechanicalVariables() 
               {
                  mpProperties=NULL;
                  mpGeometry=NULL;

                  mpB = NULL;
                  mpF0 = NULL;
                  mpDN_DX = NULL;
                  mpN = NULL;
                  mpVolumeForce = NULL;

                  number_of_variables = 0;
               };

               HydroMechanicalVariables( const GeometryType & rElementGeometry, const PropertiesType & rMaterialProperties) : mpGeometry(&rElementGeometry) , mpProperties(&rMaterialProperties)
               {
                  mpB = NULL;
                  mpF0 = NULL;
                  mpDN_DX = NULL;
                  mpN = NULL;
                  mpVolumeForce = NULL;
                  number_of_variables = 0;
               };

               // set pointers
               void SetBMatrix                   (const MatrixType& rBMatrix) {mpB=&rBMatrix;};
               void SetDeformationGradient       (const MatrixType& rF0)      {mpF0=&rF0;};
               void SetShapeFunctionsDerivatives (const MatrixType& rDN_DX)   {mpDN_DX=&rDN_DX;};
               void SetShapeFunctions            (const VectorType& rN)       {mpN = & rN; };
               void SetVolumeForce               (const VectorType& rVolumeForce){mpVolumeForce = & rVolumeForce; };

               // get const reference
               const GeometryType& GetGeometry()                {return *mpGeometry; };
               const PropertiesType& GetProperties()            {return *mpProperties; };
               const MatrixType& GetBMatrix()                   {return *mpB;};
               const MatrixType& GetDeformationGradient()       {return *mpF0;};
               const MatrixType& GetShapeFunctionsDerivatives() {return *mpDN_DX;};
               const VectorType& GetShapeFunctions()            {return *mpN; };
               const VectorType& GetVolumeForce()               {return *mpVolumeForce; };

         };

      public:

         WaterPressureUtilities();

         virtual ~WaterPressureUtilities() {};


         VectorType& CalculateAndAddHydromechanicalRHS( HydroMechanicalVariables & rVariables, VectorType & rRightHandSide, const VectorType & rBaseClassRHS, const double & rIntegrationWeight);
         
         VectorType & CalculateAndAddStabilization( HydroMechanicalVariables & rVariables, Vector & rRightHandSide, const double & rIntegrationWeight);

         MatrixType& CalculateAndAddHydromechanicalLHS( HydroMechanicalVariables & rVariables, MatrixType & rLeftHandSide, const MatrixType & rBaseClassLHS, const double & rIntegrationWeight);

         MatrixType & CalculateAndAddStabilizationLHS( HydroMechanicalVariables & rVariables, Matrix & rLeftHandSide, const double & rIntegrationWeight);

         void GetPermeabilityTensor( const PropertiesType & rProperties, const Matrix & rTotalF, Matrix & rK,  const double & rInitial_porosity, const unsigned int & rDimension, const double & rVolume );

      protected:


         // Get Properties 

         void GetScalingConstant( double& rScalingConstant, const PropertiesType& pProperties);

         void GetPermeabilityTensor( const PropertiesType & rProperties, const Matrix & rTotalF, Matrix & rK , const double & rInitial_porosity, const double & rVolume);


         virtual void GetVoigtSize( const unsigned int & dimension, unsigned int & voigtsize, unsigned int & principal_dimension); 

         double & CalculateStabilizationFactor( HydroMechanicalVariables & rVariables, double & rAlphaStabilization);

         virtual double CalculateVolumeChange( const GeometryType & rGeometry, const Vector & rN, const Matrix & rTotalF);

         // CALCULATE RHS 
         VectorType & CalculateMassBalance_WaterPressurePart( HydroMechanicalVariables & rVariables, VectorType & rLocalRHS, const double & rIntegrationWeight);

         virtual VectorType& CalculateMassBalance_AddDisplacementPart( HydroMechanicalVariables & rVariables, VectorType & rLocalRHS, const double & rIntegrationWeight);

         VectorType& CalculateWaterInternalForcesContribution( HydroMechanicalVariables & rVariables, VectorType& rRightHandSideVector, const double & rIntegrationWeight);

         VectorType & CalculateVolumeForcesContribution( HydroMechanicalVariables & rVariables, VectorType & rRightHandSideVector, const double & rIntegrationWeight);

         VectorType & CalculateStabilizationRHS( HydroMechanicalVariables & rVariables, VectorType & rRightHandSideVector, const double & rIntegrationWeight);

         // RESHAPCE RHS
         VectorType& AddReshapeBaseClassRHS( VectorType & rRightHandSideVector, const VectorType& rBaseClassRHS, const unsigned int & number_of_variables, const unsigned int & number_of_nodes);

         VectorType& AddReshapeWaterPressureForces( VectorType & rRightHandSide, const VectorType& rPartialRHS, const unsigned int number_of_variables, const unsigned int number_of_points);

         VectorType& AddReshapeWaterInternalForcesContribution( VectorType & rRightHandSideVector, const VectorType& rPartialRHS, const unsigned int number_of_variables, const unsigned int number_of_nodes, const unsigned int dimension);


         // CALCULATE LHS
         MatrixType & ComputeWaterPressureKuug( HydroMechanicalVariables & rVariables, MatrixType & rLocalLHS, const double & rIntegrationWeight);

         virtual MatrixType & ComputeWaterPressureKUwP( HydroMechanicalVariables & rVariables, MatrixType & rLocalLHS, const double & rIntegrationWeight);

         MatrixType & ComputeWaterPressureKwPwP( HydroMechanicalVariables & rVariables, MatrixType & rLocalLHS, const double & rIntegrationWeight);

         virtual MatrixType & ComputeSolidSkeletonDeformationMatrix(HydroMechanicalVariables & rVariables, MatrixType & rLocalLHS, const double & rIntegrationWeight);

         MatrixType & ComputeDarcyFlowGeometricTerms(HydroMechanicalVariables & rVariables, MatrixType & rLocalLHS, const double & rIntegrationWeight);

         virtual MatrixType & ComputeDensityChangeTerm( HydroMechanicalVariables & rVariables, MatrixType & rLocalLHS, const double & rIntegrationWeight);

         MatrixType & CalculateStabilizationLHS( HydroMechanicalVariables & rVariables, MatrixType & rLocalLHS, const double & rIntegrationWeight);


         // RESHAPE LHS
         MatrixType & AddReshapeBaseClassLHS( MatrixType & rLeftHandSide, const MatrixType & rBaseClassLHS, const unsigned int dimension, const unsigned int number_of_variables, const unsigned int number_of_nodes);

         virtual MatrixType & AddReshapeSolidSkeletonDeformationMatrix( MatrixType & rLeftHandSide, const MatrixType & rBaseClassLHS, const unsigned int dimension, const unsigned int number_of_variables, const unsigned int number_of_nodes);

         MatrixType & AddReshapeKUU( MatrixType & rLeftHandSide, const MatrixType & rBaseClassLHS, const unsigned int dimension, const unsigned int number_of_variables, const unsigned int number_of_nodes);

         MatrixType & AddReshapeKUwP( MatrixType & rLeftHandSide, const MatrixType & rBaseClassLHS, const unsigned int dimension, const unsigned int number_of_variables, const unsigned int number_of_nodes);

         MatrixType & AddReshapeKwPwP( MatrixType & rLeftHandSide, const MatrixType & rBaseClassLHS, const unsigned int dimension, const unsigned int number_of_variables, const unsigned int number_of_nodes);

         MatrixType & AddReshapeKwPU( MatrixType & rLeftHandSide, const MatrixType & rBaseClassLHS, const unsigned int dimension, const unsigned int number_of_variables, const unsigned int number_of_nodes);


         // VARIABLES
         // TO CHOOSE THE STABILIZATION 
         bool mPPP;
         
   }; // end Class WaterPressureUtilities

} // end namespace kratos

#endif // KRATOS_WATER_PRESSURE_UTILITIES
