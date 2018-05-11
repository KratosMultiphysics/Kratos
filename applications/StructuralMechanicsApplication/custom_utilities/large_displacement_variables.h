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

#if !defined(KRATOS_LARGE_DISPLACEMENT_VARIABLES_H_INCLUDED )
#define  KRATOS_LARGE_DISPLACEMENT_VARIABLES_H_INCLUDED

// System includes

// External include

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "utilities/matrix_vector_utilities.h"

namespace Kratos
{

/**
 * @class LargeDisplacementDifferentialVariables
 * @ingroup StructuralMechanicsApplication
 * @brief Differential variables for large displacements with reference configuration.
 * @details Provides deformation gradient, shape functions gradients, and
 * Jacobian information for a geometry relative to its reference configuration.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LargeDisplacementDifferentialVariables
{
public:
    typedef Geometry<Node<3>> GeometryType;

    explicit LargeDisplacementDifferentialVariables(GeometryType const& rGeom);

    /// Deformation gradient.
    const Matrix& F(std::size_t IntegrationIndex, bool IsAxisymmetric = false);

    double DetF(std::size_t IntegrationIndex, bool IsAxisymmetric = false);

    /// Shape functions gradients in reference configuration.
    const Matrix& DN_DX0(std::size_t IntegrationIndex);

    double DetJ0(std::size_t IntegrationIndex);

    const GeometryType& GetGeometry() const
    {
        return mrGeom;
    }

private:
    void RecalculateInternals(std::size_t IntegrationIndex);

    void CalculateF(std::size_t IntegrationIndex, bool IsAxisymmetric = false);

    void CalculateDetF(std::size_t IntegrationIndex, bool IsAxisymmetric = false);

    void CalculateDN_DX0(std::size_t IntegrationIndex);

    void CalculateAxisymmetricF(std::size_t IntegrationIndex);

    GeometryType const& mrGeom;
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

/**
 * @class LargeDisplacementKinematicVariables
 * @ingroup StructuralMechanicsApplication
 * @brief Kinematic variables for large displacements with reference configuration.
 * @details Provides strain-related fields.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LargeDisplacementKinematicVariables
{
public:    
    LargeDisplacementKinematicVariables(LargeDisplacementDifferentialVariables& rDiffVars,
                                        bool IsAxisymmetric = false);

    const Matrix& StrainTensor(std::size_t IntegrationIndex);

    /* Return mutable vector for compatibility with ConstitutiveLaw. */
    Vector& StrainVector(std::size_t IntegrationIndex);

    const Matrix& B(std::size_t IntegrationIndex);

private:
    LargeDisplacementDifferentialVariables& mrDiffVars;
    bool mIsAxisymmetric;
    Matrix mStrainTensor;
    Vector mStrainVector;
    Matrix mB;
    std::size_t mCurrentIntegrationIndex = -1;
    std::size_t mStrainVectorIntegrationIndex = -1;

    void Recalculate(std::size_t IntegrationIndex);

    void CalculateStrainTensor(std::size_t IntegrationIndex);

    void CalculateB(std::size_t IntegrationIndex);

    void CalculateB_Axisymmetric(std::size_t IntegrationIndex);
};

/**
 * @class LargeDisplacementKinematicVariables
 * @ingroup StructuralMechanicsApplication
 * @brief Variables related to large displacement deformations.
 * @details @see LargeDisplacementDifferentialVariables, @see LargeDisplacementKinematicVariables.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LargeDisplacementDeformationVariables
{
public:
    typedef Geometry<Node<3>> GeometryType;

    LargeDisplacementDeformationVariables(GeometryType const& rGeom,
                                          bool IsAxisymmetric = false);

    const Matrix& F(std::size_t IntegrationIndex);

    double DetF(std::size_t IntegrationIndex);

    const Matrix& DN_DX0(std::size_t IntegrationIndex);

    double DetJ0(std::size_t IntegrationIndex);

    const Matrix& StrainTensor(std::size_t IntegrationIndex);

    Vector& StrainVector(std::size_t IntegrationIndex);

    const Matrix& B(std::size_t IntegrationIndex);

private:
    LargeDisplacementDifferentialVariables mDiffVars;
    LargeDisplacementKinematicVariables mKinVars;
    const bool mIsAxisymmetric;
};

void CalculateB_2D(Matrix const& rF, Matrix const& rDN_DX0, Matrix& rB);

void CalculateB_3D(Matrix const& rF, Matrix const& rDN_DX0, Matrix& rB);

} // namespace Kratos.
#endif // KRATOS_LARGE_DISPLACEMENT_VARIABLES_H_INCLUDED  defined 
