//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.1 $
//
//

#if !defined(KRATOS_DEVIATORIC_SHAPE_UTILITIES)
#define KRATOS_DEVIATORIC_SHAPE_UTILITIES
#define PI 3.1415926535898

/*#ifdef FIND_MAX
#undef FIND_MAX
#endif

#define FIND_MAX(a, b) ((a)>(b)?(a):(b))
*/


// System includes
#include <cmath>
#include <set>


// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{


   class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) DeviatoricPlaneShapeUtilities
   {

      public:

         typedef Matrix MatrixType;

         typedef Vector VectorType;

         typedef unsigned int IndexType;

         typedef unsigned int SizeType;

         void CalculateThirdInvariantShape( const double & rLodeAngle, const double& rFriction, const double& rSmoothingAngle, double & rEffect ) 
         {

            if ( fabs( rLodeAngle) < rSmoothingAngle)
            {
               rEffect = std::cos ( rLodeAngle) - 1.0/sqrt(3.0) * std::sin( rLodeAngle) * std::sin( rFrictionAngle) ; 
            }
            else
            {
               double A, B, C; 
               ComputeSmoothingConstants( rLodeAngle, rSmoothingAngle, rFrictionAngle, A, B, C);
               double sin = std::sin( 3.0 * rLodeAngle); 

               rEffect = A + B * sin + C * sin * sin; 

            }


         }
         void CalculateThirdInvariantShapeAddim( const double & rLodeAngle, const double& rFriction, const double& rSmoothingAngle, double & rEffect ) 
         {
            CalculateThirdInvariantShapeAddim( rLodeAngle, rFriction, rSmoothingAngle, rEffect);
            rEffect /= sqrt(3.0) /6.0 * (3.0 - std::sin( rFriction) ); 
         }

      protected:

         ComputeSmoothingConstants( rLodeAngle, rSmoothingAngle,  rFrictionAngle, A, B, C)
         {
            double sign = 1.0; 
            if ( rLodeAngle < 0)
               sign = 1.0;

            double tan = std::tan(rSmoothingAngle);
            double tan3 = std::tan(3.0rSmoothingAngle);
            double sin = std::sin(rSmoothingAngle);
            double sin6 = std::sin(6.0rSmoothingAngle);
            double cos = std::cos(rSmoothingAngle);
            double cos3 = std::cos(3.0*rSmoothingAngle);
            double cos6 = std::cos(6.0*rSmoothingAngle);

            double denom = 18.0 * pow(cos3, 3.0);

            double B1, B2;

            B1 = cos * sin6 - 6.0 * cos6 * sin;
            B2 = - ( sin * sin6 + 6 * cos6 * cos);

            B1 /= denom; 
            B2 /= denom * sqrt(3.0);

            B = B1 * sign + B2 * std::sin( rFrictionAngle);

            double C1 = - cos3*cos - 3.0*sin3 * sin; 
            double C2 = cos3 * sin - 3.0*sin3 * cos; 

            C1 /= denom;
            C2 /= denom * sqrt(3.0);
            C = C1 + C2 * sign * std::sin(rFrictionAngle);


            double A1, A2; 
            A1 = cos - B1*sin3 - C1 *sin3 * sin3; 
            A2 = -1.0/sqrt(3.0) * sin - B2 * sin3 + C2 * sin3 * sin3; 

            A = A1 + A2 *sign * std::( rFrictionAngle);

         }



   }; // end Class DeviatoricPlaneShapeUtilities

} // end ramespace Kratos

#endif // KRATOS_DEVIATORIC_SHAPE

