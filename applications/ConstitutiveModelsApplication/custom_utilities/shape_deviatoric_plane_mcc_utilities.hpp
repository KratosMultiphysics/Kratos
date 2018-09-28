//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                 LlMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_DEVIATORIC_SHAPE_UTILITIES)
#define KRATOS_DEVIATORIC_SHAPE_UTILITIES


// System includes
#include <cmath>
#include <set>


// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/variables.h"

#include "constitutive_model_utilities.hpp"

namespace Kratos
{


   class ShapeAtDeviatoricPlaneMCCUtility
   {

      public:

         typedef BoundedMatrix<double,3,3> MatrixType;

         typedef array_1d<double, 6> VectorType;

         typedef unsigned int IndexType;

         typedef unsigned int SizeType;

         static inline double & EvaluateEffectDueToThirdInvariant( double & rEffect, const double & rLodeAngle, const double & rFriction)
         {
            KRATOS_TRY

            rEffect = 1.0;
            if ( rFriction < 1e-6)
               return rEffect;

            double Friction = rFriction * Globals::Pi / 180.0;
            double LodeCut = GetSmoothingLodeAngle( );
            if ( fabs( rLodeAngle) < LodeCut)
            {
               rEffect = std::cos( rLodeAngle) - 1.0/sqrt(3.0) * std::sin(Friction) * std::sin(rLodeAngle);
            } else {
               double A, B;
               GetSmoothingConstants(A, B, rLodeAngle, Friction);
               rEffect = A + B*std::sin(3.0*rLodeAngle);
            }
            rEffect /= ( sqrt(3)/6) * (3.0 - std::sin(Friction) );
            return rEffect;

            KRATOS_CATCH("")

         }

      protected:

         static inline void GetSmoothingConstants(double& rA, double& rB, const double& rLodeAngle, const double & rFriction)
         {


            double SmoothingAngle = GetSmoothingLodeAngle();

            double Sign = 1.0;
            if ( rLodeAngle < 0.0)
               Sign = -1.0;

            rA = 3.0 +  std::tan(SmoothingAngle) * std::tan(3.0*SmoothingAngle) + Sign * (std::tan( 3.0*SmoothingAngle) - 3.0*std::tan(SmoothingAngle)) * std::sin( rFriction) / sqrt(3.0);
            rA *= (1.0/3.0) * std::cos( SmoothingAngle );

            rB = -1.0 * ( Sign* std::sin(SmoothingAngle) + std::sin(rFriction)*std::cos(SmoothingAngle) / sqrt(3.0) ) / ( 3.0*std::cos(3.0*SmoothingAngle) );



         }

         static inline double GetSmoothingLodeAngle()
         {
            return 27.0*Globals::Pi/180.0;
         }



   }; // end Class ShapeAtDeviatoricPlaneMCCUtility

} // end namespace Kratos

#endif // KRATOS_DEVIATORIC_SHAPE_UTILITIES
