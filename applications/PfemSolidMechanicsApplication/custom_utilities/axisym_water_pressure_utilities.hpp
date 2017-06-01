//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.1 $
//
//

#if !defined(KRATOS_AXISYM_WATER_PRESSURE_UTILITIES)
#define KRATOS_AXISYM_WATER_PRESSURE_UTILITIES


// System includes

// External includes
//#include "boost/smart_ptr.hpp"

// Project includes
#include "custom_utilities/water_pressure_utilities.hpp" 

namespace Kratos
{


   class AxisymWaterPressureUtilities
      : public WaterPressureUtilities
   {

      public:
         typedef Matrix MatrixType;

         typedef Vector VectorType;

         typedef unsigned int IndexType;

         typedef unsigned int SizeType;

         typedef Properties PropertiesType; 

         typedef Node < 3 > NodeType; 
         typedef Geometry <NodeType>   GeometryType; 

         AxisymWaterPressureUtilities();

         virtual ~AxisymWaterPressureUtilities() {};




      protected:

         // Get Properties 
         virtual void GetVoigtSize( const unsigned int dimension, unsigned int & voigtsize, unsigned int & principal_dimension); 

         // COMPUTE RHS
         virtual VectorType& CalculateAndAddWaterPressureForcesDisplacement( VectorType& rRightHandSide , GeometryType & rGeometry,  const PropertiesType & rProperties, const MatrixType & rDN_DX, const Vector & rN, const double & rDetF0, const Matrix & rTotalF, const double & rDeltaTime, const double & rIntegrationWeight, const double & rCurrentRadius);
         
        // COMPUTE LHS


         virtual MatrixType & ComputeWaterPressureKUwP( MatrixType & LocalLHS, GeometryType & rGeometry, const Matrix & rDN_DX, const VectorType & rN, const double & rIntegrationWeight, const double & rCurrentRadius);

         virtual MatrixType & ComputeSolidSkeletonDeformationMatrix(MatrixType & rLocalLHS, GeometryType & rGeometry, const PropertiesType & rProperties, const Matrix & rDN_DX, const Vector & rN, const double & rIntegrationWeight, const double & rCurrentRadius);


         virtual MatrixType & ComputeDensityChangeTerm( MatrixType & rLocalLHS, GeometryType & rGeometry, const PropertiesType & rProperties, const Vector & rVolumeForce,  const Matrix & rDN_DX, const Vector & rN, const double & rDetF0, const double & rIntegrationWeight, const double & rCurrentRadius);

   }; // end Class AxisymWaterPressureUtilities

} // end namespace kratos

#endif // KRATOS_AXISYM_WATER_PRESSURE_UTILITIES
