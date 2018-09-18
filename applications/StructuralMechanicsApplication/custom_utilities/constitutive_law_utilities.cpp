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
    const array_1d<double, VoigtSize>& rStressVector,
    double& rI1
    )
{
    rI1 = rStressVector[0];
    for (IndexType i = 1; i < Dimension; ++i) 
        rI1 += rStressVector[i];
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateI2Invariant(
    const array_1d<double, VoigtSize>& rStressVector,
    double& rI2
    )
{
    rI2 = (rStressVector[0] + rStressVector[2]) * rStressVector[1] + rStressVector[0] * rStressVector[2] +
            -rStressVector[3] * rStressVector[3] - rStressVector[4] * rStressVector[4] - rStressVector[5] * rStressVector[5];
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateI3Invariant(
    const array_1d<double, VoigtSize>& rStressVector,
    double& rI3
    )
{
    rI3 = (rStressVector[1] * rStressVector[2] - rStressVector[4] * rStressVector[4]) * rStressVector[0] -
            rStressVector[1] * rStressVector[5] * rStressVector[5] - rStressVector[2] * rStressVector[3] * rStressVector[3] +
            2.0 * rStressVector[3] * rStressVector[4] * rStressVector[5];
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateJ2Invariant(
    const array_1d<double, VoigtSize>& rStressVector,
    const double I1,
    array_1d<double, VoigtSize>& rDeviator,
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
    for (IndexType i = Dimension; i < TVoigtSize; ++i)
        rJ2 += std::pow(rDeviator[i], 2);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateJ3Invariant(
    const array_1d<double, VoigtSize>& rDeviator,
    double& rJ3
    )
{
    rJ3 = rDeviator[0] * (rDeviator[1] * rDeviator[2] - rDeviator[4] * rDeviator[4]) +
            rDeviator[3] * (-rDeviator[3] * rDeviator[2] + rDeviator[5] * rDeviator[4]) +
            rDeviator[5] * (rDeviator[3] * rDeviator[4] - rDeviator[5] * rDeviator[1]);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateFirstVector(array_1d<double, VoigtSize>& rFirstVector)
{
    rFirstVector = ZeroVector(TVoigtSize);
    for (IndexType i = 0; i < Dimension; ++i) 
        rFirstVector[i] = 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateSecondVector(
    const array_1d<double, VoigtSize>& rDeviator,
    const double J2,
    array_1d<double, VoigtSize>& rSecondVector
    )
{
    if (rSecondVector.size() != TVoigtSize)
        rSecondVector.resize(TVoigtSize);
    const double twosqrtJ2 = 2.0 * std::sqrt(J2);
    for (IndexType i = 0; i < TVoigtSize; ++i) {
        rSecondVector[i] = rDeviator[i] / (twosqrtJ2);
    }

    for (IndexType i = Dimension; i < TVoigtSize; ++i)
        rSecondVector[i] *= 2.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateThirdVector(
    const array_1d<double, VoigtSize>& rDeviator,
    const double J2,
    array_1d<double, VoigtSize>& rThirdVector
    )
{
    if (rThirdVector.size() != TVoigtSize)
        rThirdVector.resize(TVoigtSize);

    const double J2thirds = J2 / 3.0; // static_cast<double>(Dimension);

    rThirdVector[0] = rDeviator[1] * rDeviator[2] - rDeviator[4] * rDeviator[4] + J2thirds;
    rThirdVector[1] = rDeviator[0] * rDeviator[2] - rDeviator[5] * rDeviator[5] + J2thirds;
    rThirdVector[2] = rDeviator[0] * rDeviator[1] - rDeviator[3] * rDeviator[3] + J2thirds;
    rThirdVector[3] = 2.0 * (rDeviator[4] * rDeviator[5] - rDeviator[3] * rDeviator[2]);
    rThirdVector[4] = 2.0 * (rDeviator[3] * rDeviator[4] - rDeviator[1] * rDeviator[5]);
    rThirdVector[5] = 2.0 * (rDeviator[5] * rDeviator[3] - rDeviator[0] * rDeviator[4]);
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
    double sint3 = (-3.0 * std::sqrt(3.0) * J3) / (2.0 * J2 * std::sqrt(J2));
    if (sint3 < -0.95)
        sint3 = -1.0;
    if (sint3 > 0.95)
        sint3 = 1.0;
    rLodeAngle = std::asin(sint3) / 3.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculatePrincipalStresses(
    array_1d<double, Dimension>& rPrincipalStressVector,
    const array_1d<double, VoigtSize>& rStressVector
    )
{
    double I1, I2, I3;
    CalculateI1Invariant(rStressVector, I1);
    CalculateI2Invariant(rStressVector, I2);
    CalculateI3Invariant(rStressVector, I3);
    const double II1 = std::pow(I1, 2);

    const double R = (2.0 * II1 - 9.0 * I2 * I1 + 27.0 * I3)/54.0;
    const double Q = (3.0 * I2 - II1)/9.0;

    if (std::abs(Q) > tolerance) {
        const double phi = std::acos(R / (std::sqrt(-std::pow(Q, 3))));
        const double phi_3 = phi/3.0;

        const double aux1 = 2.0 * std::sqrt(-Q);
        const double aux2 = I1 / 3.0;
        const double deg_120 = 2.0/3.0 * Globals::Pi;
        const double deg_240 = 2 * deg_120;

        rPrincipalStressVector[0] = aux2 + aux1 * std::cos(phi_3);
        rPrincipalStressVector[1] = aux2 + aux1 * std::cos(phi_3 + deg_120);
        rPrincipalStressVector[2] = aux2 + aux1 * std::cos(phi_3 + deg_240);
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
