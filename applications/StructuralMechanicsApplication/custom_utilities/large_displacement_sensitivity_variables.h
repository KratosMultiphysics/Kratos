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

#if !defined(KRATOS_LARGE_DISPLACEMENT_SENSITIVITY_VARIABLES_H_INCLUDED )
#define  KRATOS_LARGE_DISPLACEMENT_SENSITIVITY_VARIABLES_H_INCLUDED

// System includes

// External include

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "utilities/matrix_vector_utilities.h"
#include "custom_utilities/large_displacement_variables.h"
#include "utilities/geometrical_sensitivity_utility.h"

namespace Kratos
{

/**
 * @class LargeDisplacementDifferentialSensitivityVariables
 * @ingroup StructuralMechanicsApplication
 * @brief Analytical shape sensitivities for large displacement differential variables.
 * @details @see LargeDisplacementDifferentialVariables.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LargeDisplacementDifferentialSensitivityVariables
{
    public:
        typedef Geometry<Node<3>> GeometryType;

        explicit LargeDisplacementDifferentialSensitivityVariables(GeometryType const& rGeom);

        double DetJ0(std::size_t IntegrationIndex, ShapeParameter Deriv);

        const Matrix& DN_DX0(std::size_t IntegrationIndex, ShapeParameter Deriv);

        /// Deformation gradient sensitivity.
        const Matrix& F(std::size_t IntegrationIndex,
                        ShapeParameter Deriv,
                        bool IsAxisymmetric = false);

    private:
        GeometryType const& mrGeom;
        Matrix mJ0;
        Matrix mDN_DX0_Deriv;
        Matrix mF_Deriv;
        double mDetJ0_Deriv;
        std::size_t mCurrentIntegrationIndex = -1;
        ShapeParameter mCurrentShapeDerivative;

        void Synchronize(std::size_t IntegrationIndex, ShapeParameter Deriv);

        void RecalculateJ0(std::size_t IntegrationIndex);

        void RecalculateSensitivities(std::size_t IntegrationIndex, ShapeParameter Deriv);

        void CalculateFSensitivity();
};

/**
 * @class LargeDisplacementKinematicSensitivityVariables
 * @ingroup StructuralMechanicsApplication
 * @brief Analytical shape sensitivities for large displacement kinematic variables.
 * @details @see LargeDisplacementKinematicVariables.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LargeDisplacementKinematicSensitivityVariables
{
public:
    LargeDisplacementKinematicSensitivityVariables(
        LargeDisplacementDifferentialVariables& rDiffVars,
        LargeDisplacementDifferentialSensitivityVariables& rDiffVarSensitivities,
        bool IsAxisymmetric = false);

    const Matrix& StrainTensor(std::size_t IntegrationIndex, ShapeParameter Deriv);

    /* Return mutable vector for compatibility with ConstitutiveLaw. */
    Vector& StrainVector(std::size_t IntegrationIndex, ShapeParameter Deriv);

    const Matrix& B(std::size_t IntegrationIndex, ShapeParameter Deriv);

private:
    LargeDisplacementDifferentialVariables& mrDiffVars;
    LargeDisplacementDifferentialSensitivityVariables& mrDiffVarSensitivities;
    bool mIsAxisymmetric;
    Matrix mStrainTensorSensitivity;
    Vector mStrainVectorSensitivity;
    Matrix mBSensitivity;
    std::size_t mCurrentIntegrationIndex = -1;
    ShapeParameter mCurrentShapeDerivative;

    void Synchronize(std::size_t IntegrationIndex, ShapeParameter Deriv);

    void CalculateStrainTensorSensitivity(std::size_t IntegrationIndex, ShapeParameter Deriv);

    void CalculateBSensitivity(std::size_t IntegrationIndex, ShapeParameter Deriv);
};

/**
 * @class LargeDisplacementSensitivityVariables
 * @ingroup StructuralMechanicsApplication
 * @brief Analytical sensitivities for large displacement variables.
 * @details @see LargeDisplacementDifferentialVariables, @see LargeDisplacementKinematicVariables.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LargeDisplacementSensitivityVariables
{
public:
    typedef Geometry<Node<3>> GeometryType;
    
    LargeDisplacementSensitivityVariables(GeometryType const& rGeom,
                                          bool IsAxisymmetric = false);

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

} // namespace Kratos.
#endif // KRATOS_LARGE_DISPLACEMENT_SENSITIVITY_VARIABLES_H_INCLUDED  defined 
