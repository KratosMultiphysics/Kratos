//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                 LlMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_DEVIATORIC_SHAPE_MATSUOKA_UTILITIES)
#define KRATOS_DEVIATORIC_SHAPE_MATSUOKA_UTILITIES


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"

#include "constitutive_model_utilities.hpp"

namespace Kratos
{


   class ShapeAtDeviatoricPlaneMatsuokaUtility
   {

      public:

         typedef BoundedMatrix<double,3,3> MatrixType;

         typedef array_1d<double, 6> VectorType;

         typedef unsigned int IndexType;

         typedef unsigned int SizeType;

         static inline double& EvaluateEffectDueToThirdInvariant( double & rEffect, const double & rLodeAngle, const double & rFriction)
         {
            KRATOS_TRY

            rEffect = 1.0;
            //std::cout<<"Lode: "<<rLodeAngle<<std::endl;
            if ( rFriction < 1e-6)
               return rEffect;

            double Friction = rFriction * Globals::Pi / 180.0;
            //coefficient chi depending on friction angle
            double chi = std::sin(Friction) * ( pow( std::cos(Friction), 2.0) + 8.0 );
            chi /= pow( 4.0 - pow( std::cos(Friction), 2.0), 3.0/2.0);
            //effect for given Lode angle
            rEffect = 1.0/3.0 * std::acos( chi * std::sin( -3.0*rLodeAngle) );
            rEffect = 2.0*sqrt(3)*std::cos( rEffect );
            //effect triaxial compression Lode = 30° 
            double denom  = 1.0/3.0 * std::acos( -chi * 1.0 );
            denom = 2.0*sqrt(3)*std::cos( denom );
            //normalized effect g(Lode)/g(comp)
            rEffect = rEffect/denom; 

            std::cout<<"rEffect: "<<rEffect<<std::endl;
            return rEffect;

            KRATOS_CATCH("")

         }

         static inline double& EvaluateEffectDerivative( double & rDerivative, const double & rLodeAngle, const double & rFriction)
         {
            KRATOS_TRY

            rDerivative = 1.0;
            if ( rFriction < 1e-6)
               return rDerivative;

            rDerivative = 0;
            double Friction = rFriction * Globals::Pi / 180.0;
            //coefficient chi depending on friction angle
            double chi = std::sin(Friction) * ( pow( std::cos(Friction), 2.0) + 8.0 );
            chi /= pow( 4.0 - pow( std::cos(Friction), 2.0), 3.0/2.0);

            rDerivative = -chi*std::sin(std::acos(chi*std::sin(-3.0*rLodeAngle)))*std::cos(-3.0*rLodeAngle);
            rDerivative /= std::cos(1.0/3.0*std::acos(-chi*1.0))*std::sqrt(1-std::pow(chi*std::sin(-3.0*rLodeAngle),2.0));

            //std::cout<<"rDerivative: "<<rDerivative<<std::endl;
            return rDerivative;
/*
            double t2 = rLodeAngle*3.0;
            double t3 = std::cos(Friction);
            double t4 = t3*t3;
            double t5 = t4-4.0;
            double t6 = t5*t5;
            double t7 = 3.141592653589793*(1.0/3.0);
            double t8 = std::sin(Friction);
            double t9 = t4+8.0;
            double t12 = t5*t6;
            double t10 = 1.0/sqrt(-t12);
            double t11 = std::sin(t2);
            rDerivative = -(t8*t9*t10*std::sin(t7-std::acos(t8*t9*t10*t11)*(1.0/3.0))*std::cos(t2)*1.0/sqrt(1.0/(t5*t5*t5)*(t8*t8)*(t9*t9)*(t11*t11)+1.0))/std::cos(t7-acos(t8*t9*t10*8.939966636005579E-1)*(1.0/3.0));


            return rDerivative;
*/
            KRATOS_CATCH("")

         }


   }; // end Class ShapeAtDeviatoricPlaneMatsuokaUtility

} // end namespace Kratos

#endif // KRATOS_DEVIATORIC_SHAPE_UTILITIES
