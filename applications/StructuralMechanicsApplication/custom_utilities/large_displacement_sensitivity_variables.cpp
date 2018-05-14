// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    
//

// System includes

// External includes

// Project includes
#include "custom_utilities/large_displacement_sensitivity_variables.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{
LargeDisplacementDifferentialSensitivityVariables::LargeDisplacementDifferentialSensitivityVariables(
    GeometryType const& rGeom)
    : mrGeom(rGeom)
{
    mCurrentShapeDerivative.NodeIndex = -1;
    mCurrentShapeDerivative.Direction = -1;
}

double LargeDisplacementDifferentialSensitivityVariables::DetJ0(std::size_t IntegrationIndex,
                                                                ShapeParameter Deriv)
{
    KRATOS_TRY;
    Synchronize(IntegrationIndex, Deriv);
    return mDetJ0_Deriv;
    KRATOS_CATCH("");
}

const Matrix& LargeDisplacementDifferentialSensitivityVariables::DN_DX0(std::size_t IntegrationIndex,
                                                                        ShapeParameter Deriv)
{
    KRATOS_TRY;
    Synchronize(IntegrationIndex, Deriv);
    return mDN_DX0_Deriv;
    KRATOS_CATCH("");
}

const Matrix& LargeDisplacementDifferentialSensitivityVariables::F(
    std::size_t IntegrationIndex, ShapeParameter Deriv, bool IsAxisymmetric)
{
    KRATOS_TRY;
    KRATOS_ERROR_IF(IsAxisymmetric)
        << "Axisymmetic sensitivity not supported yet.\n";
    Synchronize(IntegrationIndex, Deriv);
    return mF_Deriv;
    KRATOS_CATCH("");
}

void LargeDisplacementDifferentialSensitivityVariables::Synchronize(std::size_t IntegrationIndex,
                                                                    ShapeParameter Deriv)
{
    KRATOS_TRY;
    if (IntegrationIndex != mCurrentIntegrationIndex)
    {
        RecalculateJ0(IntegrationIndex);
        RecalculateSensitivities(IntegrationIndex, Deriv);
    }
    else if (Deriv != mCurrentShapeDerivative)
        RecalculateSensitivities(IntegrationIndex, Deriv);
    KRATOS_CATCH("");
}

void LargeDisplacementDifferentialSensitivityVariables::RecalculateJ0(std::size_t IntegrationIndex)
{
    KRATOS_TRY;
    GeometryUtils::JacobianOnInitialConfiguration(
        mrGeom, mrGeom.IntegrationPoints()[IntegrationIndex], mJ0);
    mCurrentIntegrationIndex = IntegrationIndex;
    KRATOS_CATCH("");
}

void LargeDisplacementDifferentialSensitivityVariables::RecalculateSensitivities(
    std::size_t IntegrationIndex, ShapeParameter Deriv)
{
    KRATOS_TRY;
    const Matrix& rDN_De = mrGeom.ShapeFunctionsLocalGradients()[IntegrationIndex];
    auto sensitivity_utility = GeometricalSensitivityUtility(mJ0, rDN_De);
    sensitivity_utility.CalculateSensitivity(Deriv, mDetJ0_Deriv, mDN_DX0_Deriv);
    CalculateFSensitivity();
    mCurrentShapeDerivative = Deriv;
    KRATOS_CATCH("");
}

void LargeDisplacementDifferentialSensitivityVariables::CalculateFSensitivity()
{
    KRATOS_TRY;
    const std::size_t ws_dim = mrGeom.WorkingSpaceDimension();
    if (mF_Deriv.size1() != ws_dim || mF_Deriv.size2() != ws_dim)
        mF_Deriv.resize(ws_dim, ws_dim, false);
    noalias(mF_Deriv) = ZeroMatrix(ws_dim, ws_dim);
    for (std::size_t i = 0; i < ws_dim; ++i)
        for (std::size_t j = 0; j < ws_dim; ++j)
            for (std::size_t k = 0; k < mrGeom.PointsNumber(); ++k)
                mF_Deriv(i, j) += mrGeom[k].Coordinates()[i] * mDN_DX0_Deriv(k, j);
    KRATOS_CATCH("");
}

LargeDisplacementKinematicSensitivityVariables::LargeDisplacementKinematicSensitivityVariables(
    LargeDisplacementDifferentialVariables& rDiffVars,
    LargeDisplacementDifferentialSensitivityVariables& rDiffVarSensitivities,
    bool IsAxisymmetric)
    : mrDiffVars(rDiffVars),
      mrDiffVarSensitivities(rDiffVarSensitivities),
      mIsAxisymmetric(IsAxisymmetric)
{
    mCurrentShapeDerivative.NodeIndex = -1;
    mCurrentShapeDerivative.Direction = -1;
}

const Matrix& LargeDisplacementKinematicSensitivityVariables::StrainTensor(std::size_t IntegrationIndex,
                                                                           ShapeParameter Deriv)
{
    KRATOS_TRY;
    Synchronize(IntegrationIndex, Deriv);
    return mStrainTensorSensitivity;
    KRATOS_CATCH("");
}

Vector& LargeDisplacementKinematicSensitivityVariables::StrainVector(std::size_t IntegrationIndex,
                                                                     ShapeParameter Deriv)
{
    KRATOS_TRY;
    Synchronize(IntegrationIndex, Deriv);
    return mStrainVectorSensitivity;
    KRATOS_CATCH("");
}

const Matrix& LargeDisplacementKinematicSensitivityVariables::B(std::size_t IntegrationIndex,
                                                                ShapeParameter Deriv)
{
    KRATOS_TRY;
    Synchronize(IntegrationIndex, Deriv);
    return mBSensitivity;
    KRATOS_CATCH("");
}

void LargeDisplacementKinematicSensitivityVariables::Synchronize(std::size_t IntegrationIndex,
                                                                 ShapeParameter Deriv)
{
    KRATOS_TRY;
    if (IntegrationIndex != mCurrentIntegrationIndex || Deriv != mCurrentShapeDerivative)
    {
        CalculateStrainTensorSensitivity(IntegrationIndex, Deriv);
        mStrainVectorSensitivity =
            MathUtils<double>::StrainTensorToVector(mStrainTensorSensitivity);
        CalculateBSensitivity(IntegrationIndex, Deriv);
        mCurrentIntegrationIndex = IntegrationIndex;
        mCurrentShapeDerivative = Deriv;
    }
    KRATOS_CATCH("");
}

void LargeDisplacementKinematicSensitivityVariables::CalculateStrainTensorSensitivity(
    std::size_t IntegrationIndex, ShapeParameter Deriv)
{
    KRATOS_TRY;
    const Matrix& rF = mrDiffVars.F(IntegrationIndex, mIsAxisymmetric);
    const Matrix& rF_deriv =
        mrDiffVarSensitivities.F(IntegrationIndex, Deriv, mIsAxisymmetric);
    mStrainTensorSensitivity =
        0.5 * (prod(trans(rF_deriv), rF) + prod(trans(rF), rF_deriv));
    KRATOS_CATCH("");
}

void LargeDisplacementKinematicSensitivityVariables::CalculateBSensitivity(std::size_t IntegrationIndex,
                                                                           ShapeParameter Deriv)
{
    KRATOS_TRY;
    KRATOS_ERROR_IF(mIsAxisymmetric) << "Axisymmetric not supported yet.\n";
    const Matrix& rF = mrDiffVars.F(IntegrationIndex, mIsAxisymmetric);
    const Matrix& rDN_DX0 = mrDiffVars.DN_DX0(IntegrationIndex);
    const Matrix& rF_deriv =
        mrDiffVarSensitivities.F(IntegrationIndex, Deriv, mIsAxisymmetric);
    const Matrix& rDN_DX0_deriv = mrDiffVarSensitivities.DN_DX0(IntegrationIndex, Deriv);
    if (mrDiffVars.GetGeometry().WorkingSpaceDimension() == 2)
    {
        CalculateB_2D(rF_deriv, rDN_DX0, mBSensitivity);
        Matrix tmp;
        CalculateB_2D(rF, rDN_DX0_deriv, tmp);
        mBSensitivity = mBSensitivity + tmp;
    }
    else
    {
        CalculateB_3D(rF_deriv, rDN_DX0, mBSensitivity);
        Matrix tmp;
        CalculateB_3D(rF, rDN_DX0_deriv, tmp);
        mBSensitivity = mBSensitivity + tmp;
    }
    KRATOS_CATCH("");
}

LargeDisplacementSensitivityVariables::LargeDisplacementSensitivityVariables(
    GeometryType const& rGeom, bool IsAxisymmetric)
    : mDiffVars(rGeom),
      mDiffSensitivityVars(rGeom),
      mKinSensitivityVars(mDiffVars, mDiffSensitivityVars, IsAxisymmetric),
      mIsAxisymmetric(IsAxisymmetric)
{
}

} // Namespace Kratos
