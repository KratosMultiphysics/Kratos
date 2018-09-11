// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz 
//                   Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "includes/global_variables.h"
#include "utilities/math_utils.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{
template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateI1Invariant(
    const Vector& StressVector,
    double& rI1
    )
{
    rI1 = StressVector[0];
    for (IndexType i = 1; i < Dimension; ++i) 
        rI1 += StressVector[i];
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateI2Invariant(
    const Vector& StressVector,
    double& rI2
    )
{
    rI2 = (StressVector[0] + StressVector[2]) * StressVector[1] + StressVector[0] * StressVector[2] +
            -StressVector[3] * StressVector[3] - StressVector[4] * StressVector[4] - StressVector[5] * StressVector[5];
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateI3Invariant(
    const Vector& StressVector,
    double& rI3
    )
{
    rI3 = (StressVector[1] * StressVector[2] - StressVector[4] * StressVector[4]) * StressVector[0] -
            StressVector[1] * StressVector[5] * StressVector[5] - StressVector[2] * StressVector[3] * StressVector[3] +
            2.0 * StressVector[3] * StressVector[4] * StressVector[5];
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateJ2Invariant(
    const Vector& StressVector,
    const double I1,
    Vector& rDeviator,
    double& rJ2
    )
{
    rDeviator = StressVector;
    const double p_mean = I1 / static_cast<double>(Dimension);

    for (IndexType i = 0; i < Dimension; ++i) 
        rDeviator[i] -= p_mean;

    rJ2 = 0.0;
    for (IndexType i = 0; i < Dimension; ++i)
        rJ2 += 0.5 * std::pow(rDeviator[i], 2);
    for (IndexType i = Dimension; i < TVoigtSize; ++i)
        rJ2 += std::pow(rDeviator[i], 2);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateJ3Invariant(
    const Vector& Deviator,
    double& rJ3
    )
{
    rJ3 = Deviator[0] * (Deviator[1] * Deviator[2] - Deviator[4] * Deviator[4]) +
            Deviator[3] * (-Deviator[3] * Deviator[2] + Deviator[5] * Deviator[4]) +
            Deviator[5] * (Deviator[3] * Deviator[4] - Deviator[5] * Deviator[1]);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateFirstVector(Vector& FirstVector)
{
    FirstVector = ZeroVector(TVoigtSize);
    for (IndexType i = 0; i < Dimension; ++i) 
        FirstVector[i] = 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateSecondVector(
    const Vector& Deviator,
    const double J2,
    Vector& SecondVector
    )
{
    if (SecondVector.size() != TVoigtSize)
        SecondVector.resize(TVoigtSize);
    const double twosqrtJ2 = 2.0 * std::sqrt(J2);
    for (IndexType i = 0; i < TVoigtSize; ++i) {
        SecondVector[i] = Deviator[i] / (twosqrtJ2);
    }

    for (IndexType i = Dimension; i < TVoigtSize; ++i)
        SecondVector[i] *= 2.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateThirdVector(
    const Vector& Deviator,
    const double J2,
    Vector& ThirdVector
    )
{
    if (ThirdVector.size() != TVoigtSize)
        ThirdVector.resize(TVoigtSize);

    const double J2thirds = J2 / 3.0; // static_cast<double>(Dimension);

    ThirdVector[0] = Deviator[1] * Deviator[2] - Deviator[4] * Deviator[4] + J2thirds;
    ThirdVector[1] = Deviator[0] * Deviator[2] - Deviator[5] * Deviator[5] + J2thirds;
    ThirdVector[2] = Deviator[0] * Deviator[1] - Deviator[3] * Deviator[3] + J2thirds;
    ThirdVector[3] = 2.0 * (Deviator[4] * Deviator[5] - Deviator[3] * Deviator[2]);
    ThirdVector[4] = 2.0 * (Deviator[3] * Deviator[4] - Deviator[1] * Deviator[5]);
    ThirdVector[5] = 2.0 * (Deviator[5] * Deviator[3] - Deviator[0] * Deviator[4]);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateLodeAngle(
    const double J2,
    const double J3,
    double& LodeAngle
    )
{
    double sint3 = (-3.0 * std::sqrt(3.0) * J3) / (2.0 * J2 * std::sqrt(J2));
    if (sint3 < -0.95)
        sint3 = -1.0;
    if (sint3 > 0.95)
        sint3 = 1.0;
    LodeAngle = std::asin(sint3) / 3.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculatePrincipalStresses(
    Vector& rPrincipalStressVector,
    const Vector& rStressVector
    )
{
    if (rPrincipalStressVector.size() != Dimension)
        rPrincipalStressVector.resize(Dimension, false);

    double I1, I2, I3, phi, numerator, denominator, II1;
    CalculateI1Invariant(rStressVector, I1);
    CalculateI2Invariant(rStressVector, I2);
    CalculateI3Invariant(rStressVector, I3);
    II1 = I1 * I1;

    numerator = (2.0 * II1 - 9.0 * I2) * I1 + 27.0 * I3;
    denominator = (II1 - 3.0 * I2);

    if (std::abs(denominator) > tolerance) {
        phi = numerator / (2.0 * denominator * std::sqrt(denominator));

        if (std::abs(phi) > 1.0) {
            if (phi > 0.0)
                phi = 1.0;
            else
                phi = -1.0;
        }

        const double acosphi = std::acos(phi);
        phi = acosphi / 3.0;

        const double aux1 = 2.0 / 3.0 * std::sqrt(II1 - 3.0 * I2);
        const double aux2 = I1 / 3.0;
        const double deg_120 = 2.0/3.0 * Globals::Pi;
        const double deg_240 = 2 * deg_120;

        rPrincipalStressVector[0] = aux2 + aux1 * std::cos(phi);
        rPrincipalStressVector[1] = aux2 + aux1 * std::cos(phi - deg_120);
        rPrincipalStressVector[2] = aux2 + aux1 * std::cos(phi - deg_240);
    } else {
        for (IndexType i = 0; i < Dimension; ++i) {
            rPrincipalStressVector[i] = rStressVector[i];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

// template class ConstitutiveLawUtilities<3>; // TODO: Properly define the 2D case
template class ConstitutiveLawUtilities<6>;

} // namespace Kratos
