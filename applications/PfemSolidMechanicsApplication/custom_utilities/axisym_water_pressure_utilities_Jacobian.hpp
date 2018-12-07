//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.1 $
//
//

#if !defined(KRATOS_AXISYM_WATER_PRESSURE_JACOBIAN_UTILITIES)
#define KRATOS_AXISYM_WATER_PRESSURE_JACOBIAN_UTILITIES


// System includes

// External includes
//#include "boost/smart_ptr.hpp"

// Project includes
#include "custom_utilities/axisym_water_pressure_utilities.hpp" 


namespace Kratos
{


   class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) AxisymWaterPressureJacobianUtilities
      : public AxisymWaterPressureUtilities
   {

      public:
         typedef Matrix MatrixType;

         typedef Vector VectorType;

         typedef unsigned int IndexType;

         typedef unsigned int SizeType;

         typedef Properties PropertiesType; 

         typedef Node < 3 > NodeType; 
         typedef Geometry <NodeType>   GeometryType; 

         AxisymWaterPressureJacobianUtilities();

         virtual ~AxisymWaterPressureJacobianUtilities() {};



      protected:



         virtual VectorType & CalculateMassBalance_AddDisplacementPart( HydroMechanicalVariables & rVariables, VectorType & rLocalRHS, const double & rIntegrationWeight) override;

         virtual double CalculateVolumeChange( const GeometryType & rGeometry, const Vector & rN, const Matrix & rTotalF) override;
         // CALCULATE LHS

         virtual MatrixType & ComputeSolidSkeletonDeformationMatrix(HydroMechanicalVariables & rVariables, MatrixType & rLocalLHS, const double & rIntegrationWeight) override;


         // RESHAPE LHS

         virtual MatrixType & AddReshapeSolidSkeletonDeformationMatrix( MatrixType & rLeftHandSide, const MatrixType & rBaseClassLHS, const unsigned int dimension, const unsigned int number_of_variables, const unsigned int number_of_nodes) override;

   }; // end Class WaterPressureJacobianUtilities

} // end namespace kratos

#endif // KRATOS_WATER_PRESSURE_JACOBIAN_UTILITIES

