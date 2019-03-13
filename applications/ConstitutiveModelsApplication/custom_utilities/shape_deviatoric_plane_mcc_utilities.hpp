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


   class ShapeAtDeviatoricPlaneMCCUtility
   {

      public:

         typedef BoundedMatrix<double,3,3> MatrixType;

         typedef array_1d<double, 6> VectorType;

         typedef unsigned int IndexType;

         typedef unsigned int SizeType;

         static inline double& EvaluateEffectDueToThirdInvariant( double& rEffect, const double& rLodeAngle, const double& rFriction)
         {
            KRATOS_TRY

            rEffect = 1.0;
            std::cout<<"rLodeAngle: "<<rLodeAngle<<std::endl;
            if ( rFriction < 1e-6)
               return rEffect;

            double Derivative = 0;
            CalculateKLodeCoefficients( rEffect, Derivative, rLodeAngle, rFriction);
            std::cout<<"rEffect: "<<rEffect<<std::endl;
            return rEffect;

            KRATOS_CATCH("")

         }

         static inline double& EvaluateEffectDerivative( double& rEffectDeriv, const double& rLodeAngle, const double& rFriction)
         {
            KRATOS_TRY

            double Effect = 0;
            CalculateKLodeCoefficients( Effect, rEffectDeriv, rLodeAngle, rFriction);
            return rEffectDeriv;

            KRATOS_CATCH("")

         }

         static inline void CalculateKLodeCoefficients( double& rEffLode, double& rEffLodeDeriv, const double& rLodeAngle, const double& rFriction)
         {
            KRATOS_TRY

            // calcualte K(Lode) and d_K/d_Lode
            double LodeCut = GetSmoothingLodeAngle();
            double Friction = rFriction * Globals::Pi / 180.0;
            double KComp = ( sqrt(3.0)/6.0) * (3.0 - std::sin(Friction) );
            double KLode, KLodeDeriv;

            if ( fabs(rLodeAngle)  < LodeCut) {
               KLode = std::cos(rLodeAngle) - 1.0/std::sqrt(3.0) * std::sin(Friction) * std::sin(rLodeAngle); 
               KLodeDeriv = -std::sin(rLodeAngle) - 1.0/std::sqrt(3.0) * std::sin(Friction) * std::cos(rLodeAngle);
            }
            else {
               double A, B;
               GetSmoothingConstants(A, B, rLodeAngle, Friction);
               KLode = A + B * std::sin(3.0*rLodeAngle);
               KLodeDeriv = 3.0 * B * std::cos(3.0*rLodeAngle);
            }

            //rKLode /= ( std::sqrt(3.0)/6.0) * (3.0 - std::sin(Friction) );
            //rKLodeDeriv /= ( std::sqrt(3.0)/6.0) * (3.0 - std::sin(Friction) );

            rEffLode = KLode/KComp;
            rEffLodeDeriv = KLodeDeriv/KComp;

            //std::cout<<"KComp: "<<KComp<<std::endl;
            //std::cout<<"KLode: "<<KLode<<std::endl;

            KRATOS_CATCH("")

         }

      protected:

         static inline void GetSmoothingConstants(double& rA, double& rB, const double& rLodeAngle, const double& rFriction)
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
            return 25.0*Globals::Pi/180.0;
         }



   }; // end Class ShapeAtDeviatoricPlaneMCCUtility

} // end namespace Kratos

#endif // KRATOS_DEVIATORIC_SHAPE_UTILITIES
