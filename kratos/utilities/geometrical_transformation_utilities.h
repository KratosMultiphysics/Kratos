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

#if !defined(KRATOS_GEOMETRICAL_TRANSFORMATION_UTILITIES)
#define KRATOS_GEOMETRICAL_TRANSFORMATION_UTILITIES

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"

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
 * @namespace GeometricalTransformationUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities necessaries for transformation (translation and rotation)
 * @author Aditya Ghantasala
 * @author Vicente Mataix Ferrandiz
 */
namespace GeometricalTransformationUtilities
{
    typedef std::size_t IndexType;
    typedef Matrix MatrixType;
    typedef Vector VectorType;

    /**
     * @brief Calculate the transformation matrix which translates the given vector along mDirOfTranslation by mDistance
     * @param Modulus It is the magnitude by which the translation should happen in the direction of mDirOfTranslation.
     * @param rMatrix It is the transformation matrix which will be calculated in this function. This should be of correct size (4x4).
     * @param rDirOfTranslation Direction of translation
     */
    void KRATOS_API(KRATOS_CORE) CalculateTranslationMatrix(
        const double Modulus,
        MatrixType& rMatrix,
        const DenseVector<double>& rDirOfTranslation
        );

    /**
     * @brief Calculate the transformation matrix which rotates the given vector around mAxisOfRotationVector and mCenterOfRotation by provided Theta and stores the result in rMatrix The following code is generated from MATLAB and is adapted here.
     * @see http://paulbourke.net/geometry/rotate/
     * @param Theta It is the angle of rotation about mAxisOfRotationVector and mCenterOfRotation.
     * @param rMatrix It is the transformation matrix which will be calculated in this function. This should be of correct size (4x4).
     * @param rAxisOfRotationVector The axis rotation vector
     * @param rCenterOfRotation The center of rotation
     */
    void KRATOS_API(KRATOS_CORE) CalculateRotationMatrix(
        const double Theta,
        MatrixType& rMatrix,
        const DenseVector<double>& rAxisOfRotationVector,
        const DenseVector<double>& rCenterOfRotation
        );

}; // namespace GeometricalTransformationUtilities
}  // namespace Kratos
#endif /* KRATOS_GEOMETRICAL_TRANSFORMATION_UTILITIES defined */
