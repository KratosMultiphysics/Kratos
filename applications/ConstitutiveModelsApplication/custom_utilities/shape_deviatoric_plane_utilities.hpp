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

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"

#include "constitutive_model_utilities.hpp"

namespace Kratos
{


   class ShapeAtDeviatoricPlaneUtility
   {

      public:

         typedef BoundedMatrix<double,3,3> MatrixType;

         typedef array_1d<double, 6> VectorType;

         typedef unsigned int IndexType;

         typedef unsigned int SizeType;

         static inline double& EvaluateEffectDueToThirdInvariant( double& rEffect, const double& rLodeAngle, const double& rFriction)
         {
            KRATOS_TRY

            double Derivative = 0;
            CalculateKLodeCoefficients( rEffect, Derivative, rLodeAngle, rFriction);
            return rEffect;

            KRATOS_CATCH("")

         }

         static inline void CalculateKLodeCoefficients( double& rKLode, double& rKLodeDeriv, const double& rLodeAngle, const double & rFriction)
         {
            KRATOS_TRY

            if ( rFriction < 1e-6) {
               rKLode = 1.0;
               rKLodeDeriv = 0.0;
               return;
            }

            double Friction = rFriction * Globals::Pi / 180.0;
            // calcualte K(Lode) and d_K/d_Lode
            double LodeCut = GetSmoothingLodeAngle();
            if ( fabs(rLodeAngle)  < LodeCut) {
               rKLode = std::cos(rLodeAngle) - 1.0/std::sqrt(3.0) * std::sin(Friction) * std::sin(rLodeAngle); 
               rKLodeDeriv = -std::sin(rLodeAngle) - 1.0/std::sqrt(3.0) * std::sin(Friction) * std::cos(rLodeAngle);
            }
            else {

               double A, B;
               GetSmoothingConstants(A, B, rLodeAngle, Friction);

               rKLode = A + B * std::sin(3.0*rLodeAngle);
               rKLodeDeriv = 3.0 * B * std::cos(3.0*rLodeAngle);
            }

            rKLode /= ( sqrt(3.0)/6.0) * (3.0 - std::sin(Friction) );
            rKLodeDeriv /= ( sqrt(3.0)/6.0) * (3.0 - std::sin(Friction) );

            KRATOS_CATCH("")

         }

      protected:

         static inline void GetSmoothingConstants(double& rA, double& rB, const double& rLodeAngle, const double& rFriction)
         {
            KRATOS_TRY

            double SmoothingAngle = GetSmoothingLodeAngle();

            double Sign = 1.0;
            if ( rLodeAngle < 0.0)
               Sign = -1.0;

            rA = 3.0 +  std::tan(SmoothingAngle) * std::tan(3.0*SmoothingAngle) + Sign * (std::tan( 3.0*SmoothingAngle) - 3.0*std::tan(SmoothingAngle)) * std::sin( rFriction) / sqrt(3.0);
            rA *= (1.0/3.0) * std::cos( SmoothingAngle );

            rB = -1.0 * ( Sign* std::sin(SmoothingAngle) + std::sin(rFriction)*std::cos(SmoothingAngle) / sqrt(3.0) ) / ( 3.0*std::cos(3.0*SmoothingAngle) );

            KRATOS_CATCH("")
         }

         static inline double GetSmoothingLodeAngle()
         {
            return 25.0*Globals::Pi/180.0;
         }



   }; // end Class ShapeAtDeviatoricPlaneUtility

} // end namespace Kratos

#endif // KRATOS_DEVIATORIC_SHAPE_UTILITIES
