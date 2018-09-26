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


   class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) AxisymWaterPressureUtilities
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
         virtual void GetVoigtSize( const unsigned int & dimension, unsigned int & voigtsize, unsigned int & principal_dimension) override; 

         // COMPUTE RHS
         virtual VectorType & CalculateMassBalance_AddDisplacementPart( HydroMechanicalVariables & rVariables, VectorType & rLocalRHS, const double & rIntegrationWeight) override;
         
        // COMPUTE LHS


         virtual MatrixType & ComputeWaterPressureKUwP( HydroMechanicalVariables & rVariables, MatrixType & rLocalLHS, const double & rIntegrationWeight) override;

         virtual MatrixType & ComputeSolidSkeletonDeformationMatrix( HydroMechanicalVariables & rVariables, MatrixType & rLocalLHS, const double & rIntegrationWeight) override;

         virtual MatrixType & ComputeDensityChangeTerm( HydroMechanicalVariables & rVariables, MatrixType & rLocalLHS, const double & rIntegrationWeight) override;

   }; // end Class AxisymWaterPressureUtilities

} // end namespace kratos

#endif // KRATOS_AXISYM_WATER_PRESSURE_UTILITIES
