// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrándiz
//

// System includes

// External includes


// Project includes
#include "custom_elements/total_lagrangian.h"
#include "utilities/math_utils.h"
#include "utilities/matrix_vector_utilities.h"
#include "utilities/geometry_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"
#include "utilities/geometrical_sensitivity_utility.h"

namespace Kratos
{
namespace
{
class LargeDisplacementDifferentialVariables
{
    public:
    explicit LargeDisplacementDifferentialVariables(Element::GeometryType const& rGeom) : mrGeom(rGeom) {}

    /// Deformation gradient.
    const Matrix& F(std::size_t IntegrationIndex, bool IsAxisymmetric=false)
    {
        KRATOS_TRY;
        if (IntegrationIndex != mIntegrationIndexF)
            CalculateF(IntegrationIndex, IsAxisymmetric);
        return mF;
        KRATOS_CATCH("");
    }

    double DetF(std::size_t IntegrationIndex, bool IsAxisymmetric = false)
    {
        KRATOS_TRY;
        if (IntegrationIndex != mIntegrationIndexDetF)
            CalculateDetF(IntegrationIndex, IsAxisymmetric);
        return mDetF;
        KRATOS_CATCH("");
    }

    /// Shape functions gradients in reference configuration.
    const Matrix& DN_DX0(std::size_t IntegrationIndex)
    {
        KRATOS_TRY;
        if (IntegrationIndex != mIntegrationIndexDN_DX0)
            CalculateDN_DX0(IntegrationIndex);
        return mDN_DX0;
        KRATOS_CATCH("");
    }

    double DetJ0(std::size_t IntegrationIndex)
    {
        KRATOS_TRY;
        if (IntegrationIndex != mInternalIntegrationIndex)
            RecalculateInternals(IntegrationIndex);
        KRATOS_ERROR_IF(mDetJ0 < 0.0) << " INVERTED. DETJ0: " << mDetJ0 << " at " << mrGeom.Center() << '.' << std::endl;
        return mDetJ0;
        KRATOS_CATCH("");
    }

    const Element::GeometryType& GetGeometry() const
    {
        return mrGeom;
    }

private:
    void RecalculateInternals(std::size_t IntegrationIndex)
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

    void CalculateF(std::size_t IntegrationIndex, bool IsAxisymmetric=false)
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

    void CalculateDetF(std::size_t IntegrationIndex, bool IsAxisymmetric = false)
    {
        KRATOS_TRY;
        const Matrix& rF = F(IntegrationIndex, IsAxisymmetric);
        mDetF = MathUtils<double>::Det(rF);
        mIntegrationIndexDetF = IntegrationIndex;
        KRATOS_CATCH("");
    }

    void CalculateDN_DX0(std::size_t IntegrationIndex)
    {
        KRATOS_TRY;
        if (IntegrationIndex != mInternalIntegrationIndex)
            RecalculateInternals(IntegrationIndex);
        Matrix const& rDN_De = mrGeom.ShapeFunctionsLocalGradients()[IntegrationIndex];
        GeometryUtils::ShapeFunctionsGradients(rDN_De, mInvJ0, mDN_DX0);
        mIntegrationIndexDN_DX0 = IntegrationIndex;
        KRATOS_CATCH("");
    }

    void CalculateAxisymmetricF(std::size_t IntegrationIndex)
    {
        KRATOS_TRY;
        BoundedMatrix<double, 2, 2> F2x2 = mF;
        MatrixVectorUtils::InitializeMatrix(mF, 3, 3, false);
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

    Element::GeometryType const& mrGeom;
    Matrix mJ;
    Matrix mInvJ0;
    double mDetJ0;
    Matrix mF;
    Matrix mDN_DX0;
    double mDetF;
    std::size_t mInternalIntegrationIndex = -1;
    std::size_t mIntegrationIndexF = -1;
    std::size_t mIntegrationIndexDN_DX0 = -1;
    std::size_t mIntegrationIndexDetF = -1;
};

class LargeDisplacementDifferentialSensitivityVariables
{
    public:
        explicit LargeDisplacementDifferentialSensitivityVariables(Element::GeometryType const& rGeom)
            : mrGeom(rGeom)
        {
            mCurrentShapeDerivative.NodeIndex = -1;
            mCurrentShapeDerivative.Direction = -1;
        }

        double DetJ0(std::size_t IntegrationIndex, ShapeParameter Deriv)
        {
            KRATOS_TRY;
            Synchronize(IntegrationIndex, Deriv);
            return mDetJ0_Deriv;
            KRATOS_CATCH("");
        }

        const Matrix& DN_DX0(std::size_t IntegrationIndex, ShapeParameter Deriv)
        {
            KRATOS_TRY;
            Synchronize(IntegrationIndex, Deriv);
            return mDN_DX0_Deriv;
            KRATOS_CATCH("");
        }

        /// Deformation gradient sensitivity.
        const Matrix& F(std::size_t IntegrationIndex, ShapeParameter Deriv, bool IsAxisymmetric=false)
        {
            KRATOS_TRY;
            KRATOS_ERROR_IF(IsAxisymmetric)
                << "Axisymmetic sensitivity not supported yet.\n";
            Synchronize(IntegrationIndex, Deriv);
            return mF_Deriv;
            KRATOS_CATCH("");
        }

    private:
        Element::GeometryType const& mrGeom;
        Matrix mJ0;
        Matrix mDN_DX0_Deriv;
        Matrix mF_Deriv;
        double mDetJ0_Deriv;
        std::size_t mCurrentIntegrationIndex = -1;
        ShapeParameter mCurrentShapeDerivative;

        void Synchronize(std::size_t IntegrationIndex, ShapeParameter Deriv)
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

        void RecalculateJ0(std::size_t IntegrationIndex)
        {
            KRATOS_TRY;
            GeometryUtils::JacobianOnInitialConfiguration(
                mrGeom, mrGeom.IntegrationPoints()[IntegrationIndex], mJ0);
            mCurrentIntegrationIndex = IntegrationIndex;
            KRATOS_CATCH("");
        }

        void RecalculateSensitivities(std::size_t IntegrationIndex, ShapeParameter Deriv)
        {
            KRATOS_TRY;
            const Matrix& rDN_De = mrGeom.ShapeFunctionsLocalGradients()[IntegrationIndex];
            auto sensitivity_utility = GeometricalSensitivityUtility(mJ0, rDN_De);
            sensitivity_utility.CalculateSensitivity(Deriv, mDetJ0_Deriv, mDN_DX0_Deriv);
            CalculateFSensitivity();
            mCurrentShapeDerivative = Deriv;
            KRATOS_CATCH("");
        }

        void CalculateFSensitivity()
        {
            KRATOS_TRY;
            const std::size_t ws_dim = mrGeom.WorkingSpaceDimension();
            MatrixVectorUtils::InitializeMatrix(mF_Deriv, ws_dim, ws_dim, true);
            for (std::size_t i = 0; i < ws_dim; ++i)
                for (std::size_t j = 0; j < ws_dim; ++j)
                    for (std::size_t k = 0; k < mrGeom.PointsNumber(); ++k)
                        mF_Deriv(i, j) +=
                            mrGeom[k].Coordinates()[i] * mDN_DX0_Deriv(k, j);
            KRATOS_CATCH("");
        }
};

void CalculateB_2D(Matrix const& rF, Matrix const& rDN_DX0, Matrix& rB)
{
    KRATOS_TRY;
    const std::size_t dimension = 2;
    const std::size_t strain_size = 3;
    const std::size_t number_of_nodes = rDN_DX0.size1();
    MatrixVectorUtils::InitializeMatrix(rB, strain_size, dimension * number_of_nodes, false);
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
    MatrixVectorUtils::InitializeMatrix(rB, strain_size, dimension * number_of_nodes, false);
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

class LargeDisplacementKinematicVariables
{
public:
    LargeDisplacementKinematicVariables(LargeDisplacementDifferentialVariables& rDiffVars,
                                        bool IsAxisymmetric = false)
        : mrDiffVars(rDiffVars), mIsAxisymmetric(IsAxisymmetric)
    {
    }

    const Matrix& StrainTensor(std::size_t IntegrationIndex)
    {
        KRATOS_TRY;
        if (IntegrationIndex != mCurrentIntegrationIndex)
            Recalculate(IntegrationIndex);
        return mStrainTensor;
        KRATOS_CATCH("");
    }

    /* Return mutable vector for compatibility with ConstitutiveLaw. */
    Vector& StrainVector(std::size_t IntegrationIndex)
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

    const Matrix& B(std::size_t IntegrationIndex)
    {
        KRATOS_TRY;
        if (IntegrationIndex != mCurrentIntegrationIndex)
            Recalculate(IntegrationIndex);
        return mB;
        KRATOS_CATCH("");
    }

private:
    LargeDisplacementDifferentialVariables& mrDiffVars;
    bool mIsAxisymmetric;
    Matrix mStrainTensor;
    Vector mStrainVector;
    Matrix mB;
    std::size_t mCurrentIntegrationIndex = -1;
    std::size_t mStrainVectorIntegrationIndex = -1;

    void Recalculate(std::size_t IntegrationIndex)
    {
        KRATOS_TRY;
        CalculateStrainTensor(IntegrationIndex);
        CalculateB(IntegrationIndex);
        mCurrentIntegrationIndex = IntegrationIndex;
        KRATOS_CATCH("");
    }

    void CalculateStrainTensor(std::size_t IntegrationIndex)
    {
        KRATOS_TRY;
        const Matrix& rF = mrDiffVars.F(IntegrationIndex, mIsAxisymmetric);
        MatrixVectorUtils::InitializeMatrix(mStrainTensor, rF.size2(), rF.size2(), false);
        noalias(mStrainTensor) = 0.5 * prod(trans(rF), rF);
        for (std::size_t i = 0; i < mStrainTensor.size1(); ++i)
            mStrainTensor(i, i) -= 0.5;
        KRATOS_CATCH("");
    }

    void CalculateB(std::size_t IntegrationIndex)
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

    void CalculateB_Axisymmetric(std::size_t IntegrationIndex)
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
        MatrixVectorUtils::InitializeMatrix(mB, strain_size,
                                            dimension * number_of_nodes, false);
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
};

class LargeDisplacementKinematicSensitivityVariables
{
public:
    LargeDisplacementKinematicSensitivityVariables(
        LargeDisplacementDifferentialVariables& rDiffVars,
        LargeDisplacementDifferentialSensitivityVariables& rDiffVarSensitivities,
        bool IsAxisymmetric = false)
        : mrDiffVars(rDiffVars),
          mrDiffVarSensitivities(rDiffVarSensitivities),
          mIsAxisymmetric(IsAxisymmetric)
    {
        mCurrentShapeDerivative.NodeIndex = -1;
        mCurrentShapeDerivative.Direction = -1;
    }

    const Matrix& StrainTensor(std::size_t IntegrationIndex, ShapeParameter Deriv)
    {
        KRATOS_TRY;
        Synchronize(IntegrationIndex, Deriv);
        return mStrainTensorSensitivity;
        KRATOS_CATCH("");
    }

    /* Return mutable vector for compatibility with ConstitutiveLaw. */
    Vector& StrainVector(std::size_t IntegrationIndex, ShapeParameter Deriv)
    {
        KRATOS_TRY;
        Synchronize(IntegrationIndex, Deriv);
        return mStrainVectorSensitivity;
        KRATOS_CATCH("");
    }

    const Matrix& B(std::size_t IntegrationIndex, ShapeParameter Deriv)
    {
        KRATOS_TRY;
        Synchronize(IntegrationIndex, Deriv);
        return mBSensitivity;
        KRATOS_CATCH("");
    }

private:
    LargeDisplacementDifferentialVariables& mrDiffVars;
    LargeDisplacementDifferentialSensitivityVariables& mrDiffVarSensitivities;
    bool mIsAxisymmetric;
    Matrix mStrainTensorSensitivity;
    Vector mStrainVectorSensitivity;
    Matrix mBSensitivity;
    std::size_t mCurrentIntegrationIndex = -1;
    ShapeParameter mCurrentShapeDerivative;

    void Synchronize(std::size_t IntegrationIndex, ShapeParameter Deriv)
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

    void CalculateStrainTensorSensitivity(std::size_t IntegrationIndex, ShapeParameter Deriv)
    {
        KRATOS_TRY;
        const Matrix& rF = mrDiffVars.F(IntegrationIndex, mIsAxisymmetric);
        const Matrix& rF_deriv = mrDiffVarSensitivities.F(IntegrationIndex, Deriv, mIsAxisymmetric);
        mStrainTensorSensitivity =
            0.5 * (prod(trans(rF_deriv), rF) + prod(trans(rF), rF_deriv));
        KRATOS_CATCH("");
    }

    void CalculateBSensitivity(std::size_t IntegrationIndex, ShapeParameter Deriv)
    {
        KRATOS_TRY;
        KRATOS_ERROR_IF(mIsAxisymmetric) << "Axisymmetric not supported yet.\n";
        const Matrix& rF = mrDiffVars.F(IntegrationIndex, mIsAxisymmetric);
        const Matrix& rDN_DX0 = mrDiffVars.DN_DX0(IntegrationIndex);
        const Matrix& rF_deriv = mrDiffVarSensitivities.F(IntegrationIndex, Deriv, mIsAxisymmetric);
        const Matrix& rDN_DX0_deriv =
            mrDiffVarSensitivities.DN_DX0(IntegrationIndex, Deriv);
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
};
}
class LargeDisplacementDeformationVariables
{
public:
    LargeDisplacementDeformationVariables(Element::GeometryType const& rGeom,
                                          bool IsAxisymmetric = false)
        : mDiffVars(rGeom), mKinVars(mDiffVars, IsAxisymmetric), mIsAxisymmetric(IsAxisymmetric)
    {
    }

    const Matrix& F(std::size_t IntegrationIndex)
    {
        return mDiffVars.F(IntegrationIndex, mIsAxisymmetric);
    }

    double DetF(std::size_t IntegrationIndex)
    {
        return mDiffVars.DetF(IntegrationIndex, mIsAxisymmetric);
    }

    const Matrix& DN_DX0(std::size_t IntegrationIndex)
    {
        return mDiffVars.DN_DX0(IntegrationIndex);
    }

    double DetJ0(std::size_t IntegrationIndex)
    {
        return mDiffVars.DetJ0(IntegrationIndex);
    }

    const Matrix& StrainTensor(std::size_t IntegrationIndex)
    {
        return mKinVars.StrainTensor(IntegrationIndex);
    }

    Vector& StrainVector(std::size_t IntegrationIndex)
    {
        return mKinVars.StrainVector(IntegrationIndex);
    }

    const Matrix& B(std::size_t IntegrationIndex)
    {
        return mKinVars.B(IntegrationIndex);
    }

private:
    LargeDisplacementDifferentialVariables mDiffVars;
    LargeDisplacementKinematicVariables mKinVars;
    const bool mIsAxisymmetric;
};

class LargeDisplacementSensitivityVariables
{
public:
    LargeDisplacementSensitivityVariables(Element::GeometryType const& rGeom,
                                          bool IsAxisymmetric = false)
        : mDiffVars(rGeom),
          mDiffSensitivityVars(rGeom),
          mKinSensitivityVars(mDiffVars, mDiffSensitivityVars, IsAxisymmetric),
          mIsAxisymmetric(IsAxisymmetric)
    {
    }

    double DetJ0(std::size_t IntegrationIndex, ShapeParameter Deriv)
    {
        return mDiffSensitivityVars.DetJ0(IntegrationIndex, Deriv);
    }

    const Matrix& DN_DX0(std::size_t IntegrationIndex, ShapeParameter Deriv)
    {
        return mDiffSensitivityVars.DN_DX0(IntegrationIndex, Deriv);
    }

    const Matrix& F(std::size_t IntegrationIndex, ShapeParameter Deriv)
    {
        return mDiffSensitivityVars.F(IntegrationIndex, Deriv, mIsAxisymmetric);
    }

    const Matrix& StrainTensor(std::size_t IntegrationIndex, ShapeParameter Deriv)
    {
        return mKinSensitivityVars.StrainTensor(IntegrationIndex, Deriv);
    }

    /* Return mutable vector for compatibility with ConstitutiveLaw. */
    Vector& StrainVector(std::size_t IntegrationIndex, ShapeParameter Deriv)
    {
        return mKinSensitivityVars.StrainVector(IntegrationIndex, Deriv);
    }

    const Matrix& B(std::size_t IntegrationIndex, ShapeParameter Deriv)
    {
        return mKinSensitivityVars.B(IntegrationIndex, Deriv);
    }

private:
    LargeDisplacementDifferentialVariables mDiffVars;
    LargeDisplacementDifferentialSensitivityVariables mDiffSensitivityVars;
    LargeDisplacementKinematicSensitivityVariables mKinSensitivityVars;
    const bool mIsAxisymmetric;
};

TotalLagrangian::TotalLagrangian(IndexType NewId, GeometryType::Pointer pGeometry)
    : BaseSolidElement(NewId, pGeometry)
{
    // DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangian::TotalLagrangian( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : BaseSolidElement( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer TotalLagrangian::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_shared<TotalLagrangian>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangian::~TotalLagrangian()
{
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::CalculateAll( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag 
    )
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const unsigned int dimension = r_geom.WorkingSpaceDimension();
    const unsigned int mat_size = GetGeometry().PointsNumber() * dimension;

    if (CalculateStiffnessMatrixFlag == true)
        MatrixVectorUtils::InitializeMatrix(rLeftHandSideMatrix, mat_size, mat_size, true);

    if ( CalculateResidualVectorFlag == true )
        MatrixVectorUtils::InitializeVector(rRightHandSideVector, mat_size, true);

    const auto& integration_points = r_geom.IntegrationPoints();
    LargeDisplacementDeformationVariables deformation_vars(r_geom, IsAxisymmetric());
    ConstitutiveVariables cl_vars(GetStrainSize());
    Vector N;
    for (unsigned int g = 0; g < r_geom.IntegrationPointsNumber(); ++g)
    {
        const Vector body_force = this->GetBodyForce(integration_points, g);
        
        // Calculating weights for integration on the reference configuration
        double int_to_reference_weight = GetIntegrationWeight(integration_points, g, deformation_vars.DetJ0(g)); 
        if ( dimension == 2 && this->GetProperties().Has( THICKNESS )) 
            int_to_reference_weight *= this->GetProperties()[THICKNESS];

        CalculateStressAndConstitutiveMatrix(deformation_vars, g, cl_vars, rCurrentProcessInfo);
        if (CalculateStiffnessMatrixFlag == true)
        {
            /* Material stiffness matrix */
            this->CalculateAndAddKm(rLeftHandSideMatrix, deformation_vars.B(g), cl_vars.D, int_to_reference_weight);
            /* Geometric stiffness matrix */
            this->CalculateAndAddKg(rLeftHandSideMatrix, deformation_vars.DN_DX0(g), cl_vars.StressVector, int_to_reference_weight);
        }

        if (CalculateResidualVectorFlag == true)
        {
            N = row(r_geom.ShapeFunctionsValues(), g);
            this->CalculateAndAddResidualVector(
                rRightHandSideVector, N, deformation_vars.B(g), rCurrentProcessInfo,
                body_force, cl_vars.StressVector, int_to_reference_weight);
        }
    }

    KRATOS_CATCH( "" )
}

// This function is deprecated and is only implemented for compatibility with base solid element.
void TotalLagrangian::CalculateKinematicVariables(KinematicVariables& rThisKinematicVariables,
                                                  const unsigned int PointNumber,
                                                  const GeometryType::IntegrationMethod& rIntegrationMethod)
{
    KRATOS_TRY;
    const auto& r_geom = GetGeometry();
    KRATOS_ERROR_IF(rIntegrationMethod != r_geom.GetDefaultIntegrationMethod())
        << "Non-default integration method.\n";
    noalias(rThisKinematicVariables.N) = row(r_geom.ShapeFunctionsValues(), PointNumber);
    LargeDisplacementDeformationVariables deformation_vars(r_geom, IsAxisymmetric());
    rThisKinematicVariables.detJ0 = deformation_vars.DetJ0(PointNumber);
    noalias(rThisKinematicVariables.F) = deformation_vars.F(PointNumber);
    noalias(rThisKinematicVariables.B) = deformation_vars.B(PointNumber);
    rThisKinematicVariables.detF = deformation_vars.DetF(PointNumber);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::CalculateStressAndConstitutiveMatrix(
    LargeDisplacementDeformationVariables& rDeformationVars,
    std::size_t IntegrationIndex,
    ConstitutiveVariables& rOutput,
    ProcessInfo const& rCurrentProcessInfo)
{
    KRATOS_TRY;
    ConstitutiveLaw::Parameters cl_params(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    cl_params.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS |
                               ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR |
                               ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    cl_params.SetDeformationGradientF(rDeformationVars.F(IntegrationIndex));
    cl_params.SetDeterminantF(rDeformationVars.DetF(IntegrationIndex));
    cl_params.SetStrainVector(rDeformationVars.StrainVector(IntegrationIndex));
    cl_params.SetStressVector(rOutput.StressVector);
    cl_params.SetConstitutiveMatrix(rOutput.D);
    mConstitutiveLawVector[IntegrationIndex]->CalculateMaterialResponse(
        cl_params, GetStressMeasure());
    KRATOS_CATCH("");
}

void TotalLagrangian::CalculateStress(Vector& rStrain,
                                      std::size_t IntegrationIndex,
                                      Vector& rStress,
                                      ProcessInfo const& rCurrentProcessInfo)
{
    KRATOS_TRY;
    MatrixVectorUtils::InitializeVector(rStress, rStrain.size(), false);
    ConstitutiveLaw::Parameters cl_params(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    cl_params.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS | ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    cl_params.SetStrainVector(rStrain);
    cl_params.SetStressVector(rStress);
    mConstitutiveLawVector[IntegrationIndex]->CalculateMaterialResponse(cl_params, GetStressMeasure());
    KRATOS_CATCH("");
}

std::size_t TotalLagrangian::GetStrainSize() const
{
    return GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
}

bool TotalLagrangian::IsAxisymmetric() const
{
    return (GetStrainSize() == 4);
}

void TotalLagrangian::CalculateInternalForceSensitivityContribution(
    Vector& rResidualSensitivity,
    std::size_t IntegrationIndex,
    ShapeParameter Deriv,
    LargeDisplacementDeformationVariables& rDeformationVars,
    LargeDisplacementSensitivityVariables& rSensitivityVars,
    Vector const& rStressVector,
    Vector const& rStressSensitivityVector)
{
    KRATOS_TRY;
    const auto& r_geom = GetGeometry();
    const double weight = r_geom.IntegrationPoints()[IntegrationIndex].Weight() *
                          rDeformationVars.DetJ0(IntegrationIndex);
    const double weight_deriv = r_geom.IntegrationPoints()[IntegrationIndex].Weight() *
                                rSensitivityVars.DetJ0(IntegrationIndex, Deriv);
    const Matrix& rB = rDeformationVars.B(IntegrationIndex);
    const Matrix& rB_deriv = rSensitivityVars.B(IntegrationIndex, Deriv);
    rResidualSensitivity = -weight_deriv * prod(trans(rB), rStressVector);
    rResidualSensitivity -= weight * prod(trans(rB_deriv), rStressVector);
    rResidualSensitivity -= weight * prod(trans(rB), rStressSensitivityVector);
    KRATOS_CATCH("");
}

void TotalLagrangian::CalculateAndAddExternalForceSensitivityContribution(
    Vector& rResidualSensitivity,
    std::size_t IntegrationIndex,
    ShapeParameter Deriv,
    LargeDisplacementSensitivityVariables& rSensitivityVars,
    Vector const& rN,
    Vector const& rBodyForce,
    ProcessInfo const& rCurrentProcessInfo)
{
    KRATOS_TRY;
    const auto& r_geom = GetGeometry();
    const double weight_deriv = r_geom.IntegrationPoints()[IntegrationIndex].Weight() *
                                rSensitivityVars.DetJ0(IntegrationIndex, Deriv);
    CalculateAndAddExtForceContribution(rN, rCurrentProcessInfo, rBodyForce,
                                        rResidualSensitivity, weight_deriv);
    KRATOS_CATCH("");
}
/***********************************************************************************/
/***********************************************************************************/

int TotalLagrangian::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int ier = BaseSolidElement::Check(rCurrentProcessInfo);

    return ier;

    KRATOS_CATCH( "" );
}

void TotalLagrangian::CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rDesignVariable,
                                                 Matrix& rOutput,
                                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    if (rDesignVariable == SHAPE_SENSITIVITY)
    {
        const std::size_t ws_dim = r_geom.WorkingSpaceDimension();
        const std::size_t nnodes = r_geom.PointsNumber();
        Vector residual_deriv(nnodes * ws_dim), N, body_force, stress_vector, stress_vector_deriv;
        MatrixVectorUtils::InitializeMatrix(rOutput, residual_deriv.size(),
                                            residual_deriv.size(), false);
        LargeDisplacementDeformationVariables deformation_vars(r_geom, IsAxisymmetric());
        LargeDisplacementSensitivityVariables sensitivity_vars(r_geom);
        for (std::size_t g = 0; g < r_geom.IntegrationPointsNumber(); ++g)
        {
            CalculateStress(deformation_vars.StrainVector(g), g, stress_vector,
                            rCurrentProcessInfo);
            N = row(r_geom.ShapeFunctionsValues(), g);
            body_force = GetBodyForce(r_geom.IntegrationPoints(), g);
            for (auto s = ShapeParameter::Sequence(nnodes, ws_dim); s; ++s)
            {
                const auto& deriv = s.CurrentValue();
                // Assumes constant constitutive matrix wrt design parameter.
                CalculateStress(sensitivity_vars.StrainVector(g, deriv), g,
                                stress_vector_deriv, rCurrentProcessInfo);
                CalculateInternalForceSensitivityContribution(
                    residual_deriv, g, deriv, deformation_vars,
                    sensitivity_vars, stress_vector, stress_vector_deriv);
                CalculateAndAddExternalForceSensitivityContribution(
                    residual_deriv, g, deriv, sensitivity_vars, N, body_force,
                    rCurrentProcessInfo);
                for (std::size_t k = 0; k < residual_deriv.size(); ++k)
                    rOutput(k, deriv.NodeIndex * ws_dim + deriv.Direction) =
                        residual_deriv(k);
            }
        }
    }
    else
        KRATOS_ERROR << "Unsupported variable: " << rDesignVariable << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseSolidElement );
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseSolidElement );
}

} // Namespace Kratos


