// Analytical total-strain J2 response shared by the Cut-SBM operators.
// Formulas match SolidElement::AnalyticalStress/AnalyticalConstitutiveMatrix.
#pragma once

#include <cmath>
#include <limits>
#include "includes/properties.h"
#include "includes/variables.h"

namespace Kratos
{
class AnalyticalJ2ResponseUtility
{
public:
    static void AnalyticalConstitutiveMatrix(Matrix& rD, const Vector& rStrain, const Properties& rProperties)
    {
        double nu = rProperties.GetValue(POISSON_RATIO);
        double E = rProperties.GetValue(YOUNG_MODULUS);
        // Elastic materials may carry a default zero inserted by other operators.
        const double sig_y = rProperties.Has(YIELD_STRESS) && rProperties.GetValue(YIELD_STRESS) > 0.0
            ? rProperties.GetValue(YIELD_STRESS)
            : std::numeric_limits<double>::infinity();

        double eps_xx =rStrain[0];
        double gamma_xy = rStrain[2];
        double eps_yy = rStrain[1];

        rD(0,0) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            E*(nu - 1)/((nu + 1)*(2*nu - 1))
         )
         : (
            (-E*nu*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) + (1.0/3.0)*E*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) - 4*sig_y*std::pow(2*eps_xx - eps_yy, 2)*(nu + 1)*(2*nu - 1)*std::sqrt(27*std::pow(gamma_xy, 2) + 6*std::pow(eps_xx - 2*eps_yy, 2) + 6*std::pow(eps_xx + eps_yy, 2) + 6*std::pow(2*eps_xx - eps_yy, 2)) + (4.0/3.0)*std::sqrt(3)*sig_y*(nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 3.0/2.0))/((nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2))
         ));
         rD(0,1) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            -E*nu/((nu + 1)*(2*nu - 1))
         )
         : (
            (-E*nu*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) + (1.0/3.0)*E*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) + 4*sig_y*(eps_xx - 2*eps_yy)*(2*eps_xx - eps_yy)*(nu + 1)*(2*nu - 1)*std::sqrt(27*std::pow(gamma_xy, 2) + 6*std::pow(eps_xx - 2*eps_yy, 2) + 6*std::pow(eps_xx + eps_yy, 2) + 6*std::pow(2*eps_xx - eps_yy, 2)) - 2.0/3.0*std::sqrt(3)*sig_y*(nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 3.0/2.0))/((nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2))
         ));
         rD(0,2) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            0
         )
         : (
            2*gamma_xy*sig_y*(-2*eps_xx + eps_yy)/std::pow(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2), 3.0/2.0)
         ));
         rD(1,0) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            -E*nu/((nu + 1)*(2*nu - 1))
         )
         : (
            (-E*nu*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) + (1.0/3.0)*E*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) + 4*sig_y*(eps_xx - 2*eps_yy)*(2*eps_xx - eps_yy)*(nu + 1)*(2*nu - 1)*std::sqrt(27*std::pow(gamma_xy, 2) + 6*std::pow(eps_xx - 2*eps_yy, 2) + 6*std::pow(eps_xx + eps_yy, 2) + 6*std::pow(2*eps_xx - eps_yy, 2)) - 2.0/3.0*std::sqrt(3)*sig_y*(nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 3.0/2.0))/((nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2))
         ));
         rD(1,1) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            E*(nu - 1)/((nu + 1)*(2*nu - 1))
         )
         : (
            (-E*nu*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) + (1.0/3.0)*E*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) - 4*sig_y*std::pow(eps_xx - 2*eps_yy, 2)*(nu + 1)*(2*nu - 1)*std::sqrt(27*std::pow(gamma_xy, 2) + 6*std::pow(eps_xx - 2*eps_yy, 2) + 6*std::pow(eps_xx + eps_yy, 2) + 6*std::pow(2*eps_xx - eps_yy, 2)) + (4.0/3.0)*std::sqrt(3)*sig_y*(nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 3.0/2.0))/((nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2))
         ));
         rD(1,2) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            0
         )
         : (
            2*gamma_xy*sig_y*(eps_xx - 2*eps_yy)/std::pow(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2), 3.0/2.0)
         ));
         rD(2,0) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            0
         )
         : (
            2*gamma_xy*sig_y*(-2*eps_xx + eps_yy)/std::pow(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2), 3.0/2.0)
         ));
         rD(2,1) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            0
         )
         : (
            2*gamma_xy*sig_y*(eps_xx - 2*eps_yy)/std::pow(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2), 3.0/2.0)
         ));
         rD(2,2) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            (1.0/2.0)*E/(nu + 1)
         )
         : (
            4*sig_y*(std::pow(eps_xx, 2) - eps_xx*eps_yy + std::pow(eps_yy, 2))/std::pow(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2), 3.0/2.0)
         ));
    }

    static void AnalyticalStress(Vector& rStressVector, const Vector& rStrain, const Properties& rProperties)
    {
        double nu = rProperties.GetValue(POISSON_RATIO);
        double E = rProperties.GetValue(YOUNG_MODULUS);
        // Elastic materials may carry a default zero inserted by other operators.
        const double sig_y = rProperties.Has(YIELD_STRESS) && rProperties.GetValue(YIELD_STRESS) > 0.0
            ? rProperties.GetValue(YIELD_STRESS)
            : std::numeric_limits<double>::infinity();

        double eps_xx =rStrain[0];
        double gamma_xy = rStrain[2];
        double eps_yy = rStrain[1];
        
        rStressVector[0] = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            E*(eps_xx*nu - eps_xx - eps_yy*nu)/(2*std::pow(nu, 2) + nu - 1)
         )
         : (
            (1.0/3.0)*(-E*eps_xx*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2)) - E*eps_yy*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2)) + 8*eps_xx*nu*sig_y - 4*eps_xx*sig_y - 4*eps_yy*nu*sig_y + 2*eps_yy*sig_y)/((2*nu - 1)*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2)))
         ));
         rStressVector[1] = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            E*(-eps_xx*nu + eps_yy*nu - eps_yy)/(2*std::pow(nu, 2) + nu - 1)
         )
         : (
            (1.0/3.0)*(-E*eps_xx*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2)) - E*eps_yy*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2)) - 4*eps_xx*nu*sig_y + 2*eps_xx*sig_y + 8*eps_yy*nu*sig_y - 4*eps_yy*sig_y)/((2*nu - 1)*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2)))
         ));
         rStressVector[2] = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            (1.0/2.0)*E*gamma_xy/(nu + 1)
         )
         : (
            gamma_xy*sig_y/std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))
         ));
    }
};
} // namespace Kratos
