//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//     Re-factor:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/geometrical_transformation_utilities.h"

namespace Kratos
{
namespace GeometricalTransformationUtilities
{

void CalculateTranslationMatrix(
    const double Modulus,
    MatrixType& rMatrix,
    const DenseVector<double>& rDirOfTranslation
    )
{
    noalias(rMatrix) = IdentityMatrix(4,4);
    rMatrix(0,3) = Modulus * rDirOfTranslation[0];
    rMatrix(1,3) = Modulus * rDirOfTranslation[1];
    rMatrix(2,3) = Modulus * rDirOfTranslation[2];
}

/***********************************************************************************/
/***********************************************************************************/

void CalculateRotationMatrix(
    const double Theta,
    MatrixType& rMatrix,
    const DenseVector<double>& rAxisOfRotationVector,
    const DenseVector<double>& rCenterOfRotation
    )
{
    DenseVector<double> U(3); // normalized axis of rotation
    // normalizing the axis of rotation
    double norm = 0.0;
    for (IndexType d = 0; d < 3; ++d)
        norm += rAxisOfRotationVector[d] * rAxisOfRotationVector[d];
    norm = std::sqrt(norm);
    KRATOS_ERROR_IF(norm < std::numeric_limits<double>::epsilon()) << "Norm of the provided axis of rotation is Zero !"<<std::endl;
    for (IndexType d = 0; d < 3; ++d)
        U[d] = rAxisOfRotationVector[d] / norm;

    // Constructing the transformation matrix
    const double x1 = rCenterOfRotation[0];
    const double y1 = rCenterOfRotation[1];
    const double z1 = rCenterOfRotation[2];

    const double a = U[0];
    const double b = U[1];
    const double c = U[2];

    const double t2 = std::cos(Theta);
    const double t3 = std::sin(Theta);
    const double t4 = a * a;
    const double t5 = b * b;
    const double t6 = c * c;
    const double t7 = a * b;
    const double t8 = t5 + t6;
    const double t9 = std::abs(t8) < std::numeric_limits<double>::epsilon() ? 1.0e8 : 1.0 / t8;
    const double t10 = a * c;
    const double t11 = b * t3;
    const double t12 = a * t3 * t5;
    const double t13 = a * t3 * t6;
    const double t14 = b * c * t2;
    rMatrix(0,0) = t4 + t2 * t8;
    rMatrix(0,1) = t7 - c * t3 - a * b * t2;
    rMatrix(0,2) = t10 + t11 - a * c * t2;
    rMatrix(0,3) = x1 - t4 * x1 - a * b * y1 - a * c * z1 - b * t3 * z1 + c * t3 * y1 - t2 * t5 * x1 - t2 * t6 * x1 + a * b * t2 * y1 + a * c * t2 * z1;
    rMatrix(1,0) = t7 + c * t3 - a * b * t2;
    rMatrix(1,1) = t9 * (t2 * t6 + t5 * t8 + t2 * t4 * t5);
    rMatrix(1,2) = -t9 * (t12 + t13 + t14 - b * c * t8 - b * c * t2 * t4);
    rMatrix(1,3) = -t9 * (-t8 * y1 + t2 * t6 * y1 + t5 * t8 * y1 + a * b * t8 * x1 - b * c * t2 * z1 + b * c * t8 * z1 - a * t3 * t5 * z1 - a * t3 * t6 * z1 + c * t3 * t8 * x1 + t2 * t4 * t5 * y1 - a * b * t2 * t8 * x1 + b * c * t2 * t4 * z1);
    rMatrix(2,0) = t10 - t11 - a * c * t2;
    rMatrix(2,1) = t9 * (t12 + t13 - t14 + b * c * t8 + b * c * t2 * t4);
    rMatrix(2,2) = t9 * (t2 * t5 + t6 * t8 + t2 * t4 * t6);
    rMatrix(2,3) = -t9 * (-t8 * z1 + t2 * t5 * z1 + t6 * t8 * z1 + a * c * t8 * x1 - b * c * t2 * y1 + b * c * t8 * y1 + a * t3 * t5 * y1 + a * t3 * t6 * y1 - b * t3 * t8 * x1 + t2 * t4 * t6 * z1 - a * c * t2 * t8 * x1 + b * c * t2 * t4 * y1);
    rMatrix(3,0) = 0.0;
    rMatrix(3,1) = 0.0;
    rMatrix(3,2) = 0.0;
    rMatrix(3,3) = 1.0;
}

} // namespace GeometricalTransformationUtilities
} // namespace Kratos
