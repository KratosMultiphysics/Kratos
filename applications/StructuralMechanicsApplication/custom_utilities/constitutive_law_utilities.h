// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz & Alejandro Cornejo
//

#if !defined(KRATOS_CONSTITUTIVE_LAW_UTILITIES)
#define KRATOS_CONSTITUTIVE_LAW_UTILITIES

#include "utilities/math_utils.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
/**
 * @class ConstitutiveLawUtilities
 * @ingroup StructuralMechanicsApplication
 * @brief This class includes several utilities necessaries for the computation of the constitutive law
 * @details The methods are static, so it can be called without constructing the class
 * @author Vicente Mataix Ferrandiz & Alejandro Cornejo
 */
class ConstitutiveLawUtilities
{
  public:
    ///@name Type definitions
    ///@{

    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method computes the first invariant from a given stress vector
     * @param StressVector The stress vector on Voigt notation
     * @param rI1 The first invariant
     */
    static void CalculateI1Invariant(
        const Vector &StressVector,
        double &rI1)
    {
        rI1 = StressVector[0] + StressVector[1] + StressVector[2];
    }

    /**
     * @brief This method computes the second invariant from a given stress vector
     * @param StressVector The stress vector on Voigt notation
     * @param rI2 The second invariant
     */
    static void CalculateI2Invariant(
        const Vector &StressVector,
        double &rI2)
    {
        rI2 = (StressVector[0] + StressVector[2]) * StressVector[1] + StressVector[0] * StressVector[2] +
              -StressVector[3] * StressVector[3] - StressVector[4] * StressVector[4] - StressVector[5] * StressVector[5];
    }

    /**
     * @brief This method computes the third invariant from a given stress vector
     * @param StressVector The stress vector on Voigt notation
     * @param rI3 The third invariant
     */
    static void CalculateI3Invariant(
        const Vector &StressVector,
        double &rI3)
    {
        rI3 = (StressVector[1] * StressVector[2] - StressVector[4] * StressVector[4]) * StressVector[0] -
              StressVector[1] * StressVector[5] * StressVector[5] - StressVector[2] * StressVector[3] * StressVector[3] +
              2.0 * StressVector[3] * StressVector[4] * StressVector[5];
    }

    /**
     * @brief This method computes the second invariant of J
     * @param StressVector The stress vector on Voigt notation
     * @param I1 The first invariant
     * @param Deviator The deviator of the stress
     * @param rJ2 The second invariant of J
     */
    static void CalculateJ2Invariant(
        const Vector &StressVector,
        const double I1,
        Vector &rDeviator,
        double &rJ2)
    {
        rDeviator = StressVector;
        const double Pmean = I1 / 3.0;

        rDeviator[0] -= Pmean;
        rDeviator[1] -= Pmean;
        rDeviator[2] -= Pmean;

        rJ2 = 0.5 * (rDeviator[0] * rDeviator[0] + rDeviator[1] * rDeviator[1] + rDeviator[2] * rDeviator[2]) +
              (rDeviator[3] * rDeviator[3] + rDeviator[4] * rDeviator[4] + rDeviator[5] * rDeviator[5]);
    }

    /**
     * @brief This method computes the third invariant of J
     * @param Deviator The deviator of the stress
     * @param rJ3 The third invariant of J
     */
    static void CalculateJ3Invariant(
        const Vector &Deviator,
        double &rJ3)
    {
        rJ3 = Deviator[0] * (Deviator[1] * Deviator[2] - Deviator[4] * Deviator[4]) +
              Deviator[3] * (-Deviator[3] * Deviator[2] + Deviator[5] * Deviator[4]) +
              Deviator[5] * (Deviator[3] * Deviator[4] - Deviator[5] * Deviator[1]);
    }

    /**
     * @brief This method computes the first vector
     * @param FirstVector The first vector
     */
    static void CalculateFirstVector(Vector &FirstVector)
    {
        FirstVector = ZeroVector(6);
        FirstVector[0] = 1.0;
        FirstVector[1] = 1.0;
        FirstVector[2] = 1.0;
    }

    /**
     * @brief This method computes the second vector
     * @param Deviator The deviator of the stress
     * @param J2 The resultant J2 stress
     * @param SecondVector The second vector
     */
    static void CalculateSecondVector(
        const Vector &Deviator,
        const double J2,
        Vector &SecondVector)
    {
        if (SecondVector.size() == 0)
            SecondVector.resize(6);
        const double twosqrtJ2 = 2.0 * std::sqrt(J2);
        for (int i = 0; i < 6; i++)
        {
            SecondVector[i] = Deviator[i] / (twosqrtJ2);
        }

        SecondVector[3] *= 2.0;
        SecondVector[4] *= 2.0;
        SecondVector[5] *= 2.0;
    }

    /**
     * @brief This method computes the third vector
     * @param Deviator The deviator of the stress
     * @param J2 The resultant J2 stress
     * @param ThirdVector The third vector
     */
    static void CalculateThirdVector(
        const Vector &Deviator,
        const double J2,
        Vector &ThirdVector)
    {
        ThirdVector.resize(6);
        const double J2thirds = J2 / 3.0;

        ThirdVector[0] = Deviator[1] * Deviator[2] - Deviator[4] * Deviator[4] + J2thirds;
        ThirdVector[1] = Deviator[0] * Deviator[2] - Deviator[5] * Deviator[5] + J2thirds;
        ThirdVector[2] = Deviator[0] * Deviator[1] - Deviator[3] * Deviator[3] + J2thirds;
        ThirdVector[3] = 2.0 * (Deviator[4] * Deviator[5] - Deviator[3] * Deviator[2]);
        ThirdVector[4] = 2.0 * (Deviator[3] * Deviator[4] - Deviator[1] * Deviator[5]);
        ThirdVector[5] = 2.0 * (Deviator[5] * Deviator[3] - Deviator[0] * Deviator[4]);
    }

    /**
     * @brief This method computes the lode angle
     * @param J2 The resultant J2 stress
     * @param J3 The resultant J3 stress
     * @param LodeAngle The lode angle
     */
    static void CalculateLodeAngle(
        const double J2,
        const double J3,
        double &LodeAngle)
    {
        double sint3 = (-3.0 * std::sqrt(3.0) * J3) / (2.0 * J2 * std::sqrt(J2));
        if (sint3 < -0.95)
            sint3 = -1.0;
        if (sint3 > 0.95)
            sint3 = 1.0;
        LodeAngle = std::asin(sint3) / 3.0;
    }

    /**
     * @brief This method computes the principal stresses vector
     * @param rPrincipalStressVector The vector of principal stresses
     * @param rStressVector The vector of stresses
     */
    static void CalculatePrincipalStresses(
        Vector &rPrincipalStressVector,
        const Vector &rStressVector)
    {
        if (rPrincipalStressVector.size() != 3)
            rPrincipalStressVector.resize(3, false);

        double I1, I2, I3, phi, Num, Denom, II1;
        CalculateI1Invariant(rStressVector, I1);
        CalculateI2Invariant(rStressVector, I2);
        CalculateI3Invariant(rStressVector, I3);
        II1 = I1 * I1;

        Num = (2.0 * II1 - 9.0 * I2) * I1 + 27.0 * I3;
        Denom = (II1 - 3.0 * I2);

        if (std::abs(Denom) > tolerance)
        {
            phi = Num / (2.0 * Denom * std::sqrt(Denom));

            if (std::abs(phi) > 1.0)
            {
                if (phi > 0.0)
                    phi = 1.0;
                else
                    phi = -1.0;
            }

            const double acosphi = std::acos(phi);
            phi = acosphi / 3.0;

            double aux1 = 2.0 / 3.0 * std::sqrt(II1 - 3.0 * I2);
            double aux2 = I1 / 3.0;

            rPrincipalStressVector[0] = aux2 + aux1 * std::cos(phi);
            rPrincipalStressVector[1] = aux2 + aux1 * std::cos(phi - 2.09439510239);
            rPrincipalStressVector[2] = aux2 + aux1 * std::cos(phi - 4.18879020478);
        }
        else
        {
            rPrincipalStressVector = ZeroVector(3);
        }
    }

  private:
}; // class ConstitutiveLawUtilities
} // namespace Kratos
#endif /* KRATOS_CONSTITUTIVE_LAW_UTILITIES defined */
