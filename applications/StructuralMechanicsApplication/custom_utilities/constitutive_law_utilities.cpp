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

template<>
void ConstitutiveLawUtilities<6>::CalculateI1Invariant(
    const array_1d<double, 6>& rStressVector,
    double& rI1
    )
{
    rI1 = rStressVector[0];
    for (IndexType i = 1; i < Dimension; ++i) 
        rI1 += rStressVector[i];
}

template<>
void ConstitutiveLawUtilities<3>::CalculateI1Invariant(
    const array_1d<double, 3>& rStressVector,
    double& rI1
    )
{
    rI1 = rStressVector[0];
    for (IndexType i = 1; i < Dimension; ++i) 
        rI1 += rStressVector[i];
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateI2Invariant(
    const array_1d<double, 6>& rStressVector,
    double& rI2
    )
{
    rI2 = (rStressVector[0] + rStressVector[2]) * rStressVector[1] + rStressVector[0] * rStressVector[2] +
            -rStressVector[3] * rStressVector[3] - rStressVector[4] * rStressVector[4] - rStressVector[5] * rStressVector[5];
}

template<>
void ConstitutiveLawUtilities<3>::CalculateI2Invariant(
    const array_1d<double, 3>& rStressVector,
    double& rI2
    )
{
    rI2 = rStressVector[0] * rStressVector[1] - std::pow(rStressVector[2], 2);
}
/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateI3Invariant(
    const array_1d<double, 6>& rStressVector,
    double& rI3
    )
{
    rI3 = (rStressVector[1] * rStressVector[2] - rStressVector[4] * rStressVector[4]) * rStressVector[0] -
            rStressVector[1] * rStressVector[5] * rStressVector[5] - rStressVector[2] * rStressVector[3] * rStressVector[3] +
            2.0 * rStressVector[3] * rStressVector[4] * rStressVector[5];
}

template<>
void ConstitutiveLawUtilities<3>::CalculateI3Invariant(
    const array_1d<double, 3>& rStressVector,
    double& rI3
    )
{
    KRATOS_ERROR << "I3 invariant not available in 2D!" << std::endl;
}
/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateJ2Invariant(
    const array_1d<double, 6>& rStressVector,
    const double I1,
    array_1d<double, 6>& rDeviator,
    double& rJ2
    )
{
    rDeviator = rStressVector;
    const double p_mean = I1 / static_cast<double>(Dimension);

    for (IndexType i = 0; i < Dimension; ++i) 
        rDeviator[i] -= p_mean;

    rJ2 = 0.0;
    for (IndexType i = 0; i < Dimension; ++i)
        rJ2 += 0.5 * std::pow(rDeviator[i], 2);
    for (IndexType i = Dimension; i < 6; ++i)
        rJ2 += std::pow(rDeviator[i], 2);
}

template<>
void ConstitutiveLawUtilities<3>::CalculateJ2Invariant(
    const array_1d<double, 3>& rStressVector,
    const double I1,
    array_1d<double, 3>& rDeviator,
    double& rJ2
    )
{
    rDeviator = rStressVector;
    const double p_mean = I1 / static_cast<double>(Dimension);

    for (IndexType i = 0; i < Dimension; ++i)
        rDeviator[i] -= p_mean;

    rJ2 = 0.5 * (std::pow(rDeviator[0], 2) + std::pow(rDeviator[1], 2)) +
          std::pow(rDeviator[2], 2);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateJ3Invariant(
    const array_1d<double, 6>& rDeviator,
    double& rJ3
    )
{
    rJ3 = rDeviator[0] * (rDeviator[1] * rDeviator[2] - rDeviator[4] * rDeviator[4]) +
            rDeviator[3] * (-rDeviator[3] * rDeviator[2] + rDeviator[5] * rDeviator[4]) +
            rDeviator[5] * (rDeviator[3] * rDeviator[4] - rDeviator[5] * rDeviator[1]);
}

template<>
void ConstitutiveLawUtilities<3>::CalculateJ3Invariant(
    const array_1d<double, 3>& rDeviator,
    double& rJ3
    )
{
    rJ3 = rDeviator[0] * rDeviator[1] - std::pow(rDeviator[2], 2);
}
/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateFirstVector(array_1d<double, 6>& rFirstVector)
{
    rFirstVector = ZeroVector(6);
    for (IndexType i = 0; i < Dimension; ++i) 
        rFirstVector[i] = 1.0;
}

template<>
void ConstitutiveLawUtilities<3>::CalculateFirstVector(array_1d<double, 3>& rFirstVector)
{
    // TODO
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateSecondVector(
    const array_1d<double, 6>& rDeviator,
    const double J2,
    array_1d<double, 6>& rSecondVector
    )
{
    if (rSecondVector.size() != 6)
        rSecondVector.resize(6);
    const double twosqrtJ2 = 2.0 * std::sqrt(J2);
    for (IndexType i = 0; i < 6; ++i) {
        rSecondVector[i] = rDeviator[i] / (twosqrtJ2);
    }

    for (IndexType i = Dimension; i < 6; ++i)
        rSecondVector[i] *= 2.0;
}

template<>
void ConstitutiveLawUtilities<3>::CalculateSecondVector(
    const array_1d<double, 3>& rDeviator,
    const double J2,
    array_1d<double, 3>& rSecondVector
    )
{
    // TODO
}
/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateThirdVector(
    const array_1d<double, 6>& rDeviator,
    const double J2,
    array_1d<double, 6>& rThirdVector
    )
{
    if (rThirdVector.size() != 6)
        rThirdVector.resize(6);

    const double J2thirds = J2 / 3.0; // static_cast<double>(Dimension);

    rThirdVector[0] = rDeviator[1] * rDeviator[2] - rDeviator[4] * rDeviator[4] + J2thirds;
    rThirdVector[1] = rDeviator[0] * rDeviator[2] - rDeviator[5] * rDeviator[5] + J2thirds;
    rThirdVector[2] = rDeviator[0] * rDeviator[1] - rDeviator[3] * rDeviator[3] + J2thirds;
    rThirdVector[3] = 2.0 * (rDeviator[4] * rDeviator[5] - rDeviator[3] * rDeviator[2]);
    rThirdVector[4] = 2.0 * (rDeviator[3] * rDeviator[4] - rDeviator[1] * rDeviator[5]);
    rThirdVector[5] = 2.0 * (rDeviator[5] * rDeviator[3] - rDeviator[0] * rDeviator[4]);
}

template<>
void ConstitutiveLawUtilities<3>::CalculateThirdVector(
    const array_1d<double, 3>& rDeviator,
    const double J2,
    array_1d<double, 3>& rThirdVector
    )
{
    // TODO
}
/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateLodeAngle(
    const double J2,
    const double J3,
    double& rLodeAngle
    )
{
    if (std::abs(J2) > tolerance) {
        double sint3 = (-3.0 * std::sqrt(3.0) * J3) / (2.0 * J2 * std::sqrt(J2));
        if (sint3 < -0.95)
            sint3 = -1.0;
        else if (sint3 > 0.95)
            sint3 = 1.0;
        rLodeAngle = std::asin(sint3) / 3.0;
    } else {
        rLodeAngle = 0.0;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculatePrincipalStresses(
    array_1d<double, Dimension>& rPrincipalStressVector,
    const array_1d<double, 6>& rStressVector
    )
{
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

template<>
void ConstitutiveLawUtilities<3>::CalculatePrincipalStresses(
    array_1d<double, Dimension>& rPrincipalStressVector,
    const array_1d<double, 3>& rStressVector
    )
{
    if (rPrincipalStressVector.size() != Dimension)
		rPrincipalStressVector.resize(Dimension);

    rPrincipalStressVector[0] = 0.5 * (rStressVector[0] + rStressVector[1]) + std::sqrt(std::pow(0.5 * (rStressVector[0] - rStressVector[1]), 2) + std::pow(rStressVector[2], 2));
    rPrincipalStressVector[1] = 0.5 * (rStressVector[0] + rStressVector[1]) - std::sqrt(std::pow(0.5 * (rStressVector[0] - rStressVector[1]), 2) + std::pow(rStressVector[2], 2));
}
/***********************************************************************************/
/***********************************************************************************/

template class ConstitutiveLawUtilities<3>; 
template class ConstitutiveLawUtilities<6>;







} // namespace Kratos
