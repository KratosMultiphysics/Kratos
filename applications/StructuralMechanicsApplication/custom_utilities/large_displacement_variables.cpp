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
#include "custom_utilities/large_displacement_variables.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"

namespace Kratos
{
LargeDisplacementDifferentialVariables::LargeDisplacementDifferentialVariables(GeometryType const& rGeom)
    : mrGeom(rGeom)
{
}

const Matrix& LargeDisplacementDifferentialVariables::F(std::size_t IntegrationIndex, bool IsAxisymmetric)
{
    KRATOS_TRY;
    if (IntegrationIndex != mIntegrationIndexF)
        CalculateF(IntegrationIndex, IsAxisymmetric);
    return mF;
    KRATOS_CATCH("");
}

double LargeDisplacementDifferentialVariables::DetF(std::size_t IntegrationIndex, bool IsAxisymmetric)
{
    KRATOS_TRY;
    if (IntegrationIndex != mIntegrationIndexDetF)
        CalculateDetF(IntegrationIndex, IsAxisymmetric);
    return mDetF;
    KRATOS_CATCH("");
}

const Matrix& LargeDisplacementDifferentialVariables::DN_DX0(std::size_t IntegrationIndex)
{
    KRATOS_TRY;
    if (IntegrationIndex != mIntegrationIndexDN_DX0)
        CalculateDN_DX0(IntegrationIndex);
    return mDN_DX0;
    KRATOS_CATCH("");
}

double LargeDisplacementDifferentialVariables::DetJ0(std::size_t IntegrationIndex)
{
    KRATOS_TRY;
    if (IntegrationIndex != mInternalIntegrationIndex)
        RecalculateInternals(IntegrationIndex);
    KRATOS_ERROR_IF(mDetJ0 < 0.0) << " INVERTED. DETJ0: " << mDetJ0 << " at "
                                  << mrGeom.Center() << '.' << std::endl;
    return mDetJ0;
    KRATOS_CATCH("");
}

void LargeDisplacementDifferentialVariables::RecalculateInternals(std::size_t IntegrationIndex)
{
    KRATOS_TRY;
    // Calculate mInvJ0 and mDetJ0.
    GeometryUtils::JacobianOnInitialConfiguration(
        mrGeom, mrGeom.IntegrationPoints()[IntegrationIndex], mJ);
    MathUtils<double>::InvertMatrix(mJ, mInvJ0, mDetJ0);
    // Calculate mJ.
    mrGeom.Jacobian(mJ, IntegrationIndex);
    // Update current integration point.
    mInternalIntegrationIndex = IntegrationIndex;
    KRATOS_CATCH("");
}

void LargeDisplacementDifferentialVariables::CalculateF(std::size_t IntegrationIndex,
                                                        bool IsAxisymmetric)
{
    KRATOS_TRY;
    if (IntegrationIndex != mInternalIntegrationIndex)
        RecalculateInternals(IntegrationIndex);
    GeometryUtils::DeformationGradient(mJ, mInvJ0, mF);
    if (IsAxisymmetric)
        CalculateAxisymmetricF(IntegrationIndex);
    mIntegrationIndexF = IntegrationIndex;
    KRATOS_CATCH("");
}

void LargeDisplacementDifferentialVariables::CalculateDetF(std::size_t IntegrationIndex,
                                                           bool IsAxisymmetric)
{
    KRATOS_TRY;
    const Matrix& rF = F(IntegrationIndex, IsAxisymmetric);
    mDetF = MathUtils<double>::Det(rF);
    mIntegrationIndexDetF = IntegrationIndex;
    KRATOS_CATCH("");
}

void LargeDisplacementDifferentialVariables::CalculateDN_DX0(std::size_t IntegrationIndex)
{
    KRATOS_TRY;
    if (IntegrationIndex != mInternalIntegrationIndex)
        RecalculateInternals(IntegrationIndex);
    Matrix const& rDN_De = mrGeom.ShapeFunctionsLocalGradients()[IntegrationIndex];
    GeometryUtils::ShapeFunctionsGradients(rDN_De, mInvJ0, mDN_DX0);
    mIntegrationIndexDN_DX0 = IntegrationIndex;
    KRATOS_CATCH("");
}

void LargeDisplacementDifferentialVariables::CalculateAxisymmetricF(std::size_t IntegrationIndex)
{
    KRATOS_TRY;
    BoundedMatrix<double, 2, 2> F2x2 = mF;
    if (mF.size1() != 3 || mF.size2() != 3)
        mF.resize(3, 3, false);
    for (unsigned i = 0; i < 2; ++i)
    {
        for (unsigned j = 0; j < 2; ++j)
            mF(i, j) = F2x2(i, j);
        mF(i, 2) = mF(2, i) = 0.0;
    }
    Vector N = row(mrGeom.ShapeFunctionsValues(), IntegrationIndex);
    const double current_radius =
        StructuralMechanicsMathUtilities::CalculateRadius(N, mrGeom, Current);
    const double initial_radius =
        StructuralMechanicsMathUtilities::CalculateRadius(N, mrGeom, Initial);
    mF(2, 2) = current_radius / initial_radius;
    KRATOS_CATCH("");
}

LargeDisplacementKinematicVariables::LargeDisplacementKinematicVariables(
    LargeDisplacementDifferentialVariables& rDiffVars, bool IsAxisymmetric)
    : mrDiffVars(rDiffVars), mIsAxisymmetric(IsAxisymmetric)
{
}

const Matrix& LargeDisplacementKinematicVariables::StrainTensor(std::size_t IntegrationIndex)
{
    KRATOS_TRY;
    if (IntegrationIndex != mCurrentIntegrationIndex)
        Recalculate(IntegrationIndex);
    return mStrainTensor;
    KRATOS_CATCH("");
}

Vector& LargeDisplacementKinematicVariables::StrainVector(std::size_t IntegrationIndex)
{
    KRATOS_TRY;
    if (IntegrationIndex != mStrainVectorIntegrationIndex)
    {
        mStrainVector =
            MathUtils<double>::StrainTensorToVector(StrainTensor(IntegrationIndex));
        mStrainVectorIntegrationIndex = IntegrationIndex;
    }
    return mStrainVector;
    KRATOS_CATCH("");
}

const Matrix& LargeDisplacementKinematicVariables::B(std::size_t IntegrationIndex)
{
    KRATOS_TRY;
    if (IntegrationIndex != mCurrentIntegrationIndex)
        Recalculate(IntegrationIndex);
    return mB;
    KRATOS_CATCH("");
}

void LargeDisplacementKinematicVariables::Recalculate(std::size_t IntegrationIndex)
{
    KRATOS_TRY;
    CalculateStrainTensor(IntegrationIndex);
    CalculateB(IntegrationIndex);
    mCurrentIntegrationIndex = IntegrationIndex;
    KRATOS_CATCH("");
}

void LargeDisplacementKinematicVariables::CalculateStrainTensor(std::size_t IntegrationIndex)
{
    KRATOS_TRY;
    const Matrix& rF = mrDiffVars.F(IntegrationIndex, mIsAxisymmetric);
    if (mStrainTensor.size1() != rF.size2() || mStrainTensor.size2() != rF.size2())
        mStrainTensor.resize(rF.size2(), rF.size2(), false);
    noalias(mStrainTensor) = 0.5 * prod(trans(rF), rF);
    for (std::size_t i = 0; i < mStrainTensor.size1(); ++i)
        mStrainTensor(i, i) -= 0.5;
    KRATOS_CATCH("");
}

void LargeDisplacementKinematicVariables::CalculateB(std::size_t IntegrationIndex)
{
    KRATOS_TRY;
    if (mIsAxisymmetric)
        CalculateB_Axisymmetric(IntegrationIndex);
    else if (mrDiffVars.GetGeometry().WorkingSpaceDimension() == 2)
        CalculateB_2D(mrDiffVars.F(IntegrationIndex),
                      mrDiffVars.DN_DX0(IntegrationIndex), mB);
    else
        CalculateB_3D(mrDiffVars.F(IntegrationIndex),
                      mrDiffVars.DN_DX0(IntegrationIndex), mB);
    KRATOS_CATCH("");
}

void LargeDisplacementKinematicVariables::CalculateB_Axisymmetric(std::size_t IntegrationIndex)
{
    KRATOS_TRY;
    const auto& r_geom = mrDiffVars.GetGeometry();
    const unsigned int number_of_nodes = r_geom.PointsNumber();
    const unsigned int dimension = r_geom.WorkingSpaceDimension();
    const unsigned int strain_size = 4;
    const Matrix& rF = mrDiffVars.F(IntegrationIndex, true);
    const Matrix& rDN_DX0 = mrDiffVars.DN_DX0(IntegrationIndex);
    Vector N = row(r_geom.ShapeFunctionsValues(), IntegrationIndex);
    double radius = StructuralMechanicsMathUtilities::CalculateRadius(N, r_geom);
    
    if (mB.size1() != strain_size || mB.size2() != dimension * number_of_nodes)
        mB.resize(strain_size, dimension * number_of_nodes, false);
    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        const unsigned int index = dimension * i;
        mB(0, index + 0) = rF(0, 0) * rDN_DX0(i, 0);
        mB(0, index + 1) = rF(1, 0) * rDN_DX0(i, 0);
        mB(1, index + 1) = rF(0, 1) * rDN_DX0(i, 1);
        mB(1, index + 1) = rF(1, 1) * rDN_DX0(i, 1);
        mB(2, index + 0) = N[i] / radius;
        mB(3, index + 0) = rF(0, 0) * rDN_DX0(i, 1) + rF(0, 1) * rDN_DX0(i, 0);
        mB(3, index + 1) = rF(1, 0) * rDN_DX0(i, 1) + rF(1, 1) * rDN_DX0(i, 0);
    }
    KRATOS_CATCH("");
}

LargeDisplacementDeformationVariables::LargeDisplacementDeformationVariables(
    GeometryType const& rGeom, bool IsAxisymmetric)
    : mDiffVars(rGeom), mKinVars(mDiffVars, IsAxisymmetric), mIsAxisymmetric(IsAxisymmetric)
{
}

const Matrix& LargeDisplacementDeformationVariables::F(std::size_t IntegrationIndex)
{
    return mDiffVars.F(IntegrationIndex, mIsAxisymmetric);
}

double LargeDisplacementDeformationVariables::DetF(std::size_t IntegrationIndex)
{
    return mDiffVars.DetF(IntegrationIndex, mIsAxisymmetric);
}

const Matrix& LargeDisplacementDeformationVariables::DN_DX0(std::size_t IntegrationIndex)
{
    return mDiffVars.DN_DX0(IntegrationIndex);
}

double LargeDisplacementDeformationVariables::DetJ0(std::size_t IntegrationIndex)
{
    return mDiffVars.DetJ0(IntegrationIndex);
}

const Matrix& LargeDisplacementDeformationVariables::StrainTensor(std::size_t IntegrationIndex)
{
    return mKinVars.StrainTensor(IntegrationIndex);
}

Vector& LargeDisplacementDeformationVariables::StrainVector(std::size_t IntegrationIndex)
{
    return mKinVars.StrainVector(IntegrationIndex);
}

const Matrix& LargeDisplacementDeformationVariables::B(std::size_t IntegrationIndex)
{
    return mKinVars.B(IntegrationIndex);
}

void CalculateB_2D(Matrix const& rF, Matrix const& rDN_DX0, Matrix& rB)
{
    KRATOS_TRY;
    const std::size_t dimension = 2;
    const std::size_t strain_size = 3;
    const std::size_t number_of_nodes = rDN_DX0.size1();
    if (rB.size1() != strain_size || rB.size2() != dimension * number_of_nodes)
        rB.resize(strain_size, dimension * number_of_nodes, false);
    for (std::size_t i = 0; i < number_of_nodes; ++i)
    {
        const auto index = dimension * i;
        rB(0, index + 0) = rF(0, 0) * rDN_DX0(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDN_DX0(i, 0);
        rB(1, index + 0) = rF(0, 1) * rDN_DX0(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDN_DX0(i, 1);
        rB(2, index + 0) = rF(0, 0) * rDN_DX0(i, 1) + rF(0, 1) * rDN_DX0(i, 0);
        rB(2, index + 1) = rF(1, 0) * rDN_DX0(i, 1) + rF(1, 1) * rDN_DX0(i, 0);
    }
    KRATOS_CATCH("")
}

void CalculateB_3D(Matrix const& rF, Matrix const& rDN_DX0, Matrix& rB)
{
    KRATOS_TRY;
    const std::size_t dimension = 3;
    const std::size_t strain_size = 6;
    const std::size_t number_of_nodes = rDN_DX0.size1();
    if (rB.size1() != strain_size || rB.size2() != dimension * number_of_nodes)
        rB.resize(strain_size, dimension * number_of_nodes, false);
    for (std::size_t i = 0; i < rDN_DX0.size1(); ++i)
    {
        const auto index = dimension * i;
        rB(0, index + 0) = rF(0, 0) * rDN_DX0(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDN_DX0(i, 0);
        rB(0, index + 2) = rF(2, 0) * rDN_DX0(i, 0);
        rB(1, index + 0) = rF(0, 1) * rDN_DX0(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDN_DX0(i, 1);
        rB(1, index + 2) = rF(2, 1) * rDN_DX0(i, 1);
        rB(2, index + 0) = rF(0, 2) * rDN_DX0(i, 2);
        rB(2, index + 1) = rF(1, 2) * rDN_DX0(i, 2);
        rB(2, index + 2) = rF(2, 2) * rDN_DX0(i, 2);
        rB(3, index + 0) = rF(0, 0) * rDN_DX0(i, 1) + rF(0, 1) * rDN_DX0(i, 0);
        rB(3, index + 1) = rF(1, 0) * rDN_DX0(i, 1) + rF(1, 1) * rDN_DX0(i, 0);
        rB(3, index + 2) = rF(2, 0) * rDN_DX0(i, 1) + rF(2, 1) * rDN_DX0(i, 0);
        rB(4, index + 0) = rF(0, 1) * rDN_DX0(i, 2) + rF(0, 2) * rDN_DX0(i, 1);
        rB(4, index + 1) = rF(1, 1) * rDN_DX0(i, 2) + rF(1, 2) * rDN_DX0(i, 1);
        rB(4, index + 2) = rF(2, 1) * rDN_DX0(i, 2) + rF(2, 2) * rDN_DX0(i, 1);
        rB(5, index + 0) = rF(0, 2) * rDN_DX0(i, 0) + rF(0, 0) * rDN_DX0(i, 2);
        rB(5, index + 1) = rF(1, 2) * rDN_DX0(i, 0) + rF(1, 0) * rDN_DX0(i, 2);
        rB(5, index + 2) = rF(2, 2) * rDN_DX0(i, 0) + rF(2, 0) * rDN_DX0(i, 2);
    }
    KRATOS_CATCH("");
}

} // Namespace Kratos


