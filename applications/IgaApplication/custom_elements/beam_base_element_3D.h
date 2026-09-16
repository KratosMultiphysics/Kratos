//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Max Friedrichs Dachale, Juan Ignacio Camarotti,
//                   Ricky Aristio, Alicia Knauer
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "containers/array_1d.h"

// Application includes
#include "iga_application_variables.h"

namespace Kratos
{

/**
 * @class BeamBaseElement3D
 * @brief Base class for three-dimensional IGA beam elements.
 *
 * This class provides functionality shared by three-dimensional beam
 * formulations, including degree-of-freedom handling, constitutive-law
 * initialization, beam kinematics, and rotation-related operations.
 */
class KRATOS_API(IGA_APPLICATION) BeamBaseElement3D
    : public Element
{
public:

    ///@name Type Definitions
    ///@{

    using BaseType = Element;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using SizeType = typename BaseType::SizeType;
    using IndexType = typename BaseType::IndexType;
    using Matrix3d = BoundedMatrix<double, 3, 3>;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(BeamBaseElement3D);

    ///@}
    ///@name Life Cycle
    ///@{

    BeamBaseElement3D(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    BeamBaseElement3D(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    BeamBaseElement3D()
        : Element()
    {
    }

    ~BeamBaseElement3D() override = default;

    ///@}
    ///@name Element Creation
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) const override = 0;

    Element::Pointer Create(
        IndexType NewId,
        const NodesArrayType& rThisNodes,
        PropertiesType::Pointer pProperties) const override = 0;

    ///@}
    ///@name Degrees of Freedom
    ///@{

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetValuesVector(
        Vector& rValues,
        int Step = 0) const override;

    ///@}
    ///@name Material
    ///@{

    void InitializeMaterial();

    ///@}
    ///@name Element Calculations
    ///@{

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override = 0;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override = 0;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override = 0;

    ///@}
    ///@name Input and Output
    ///@{

    int Check(
        const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}

protected:

    ///@name Protected Structures
    ///@{

    /**
     * @brief Variables used during constitutive-law calculations.
     */
    struct ConstitutiveVariables
    {
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;

        explicit ConstitutiveVariables(
            const SizeType StrainSize = 3)
        {
            StrainVector = ZeroVector(StrainSize);
            StressVector = ZeroVector(StrainSize);
            ConstitutiveMatrix = ZeroMatrix(StrainSize, StrainSize);
        }
    };

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Performs the formulation-specific element calculation.
     */
    virtual void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag) = 0;

    /**
     * @brief Computes the first variation of the beam tangent.
     */
    void ComputeTangentVariation(
        const Vector& rShapeFunctionDerivatives,
        const array_1d<double, 3>& rFirstBaseVector,
        Vector& rTangentVariation) const;

    /**
     * @brief Computes the second variation of the beam tangent.
     */
    void ComputeTangentSecondVariation(
        const Vector& rShapeFunctionDerivatives,
        const array_1d<double, 3>& rFirstBaseVector,
        Matrix& rTangentSecondVariation) const;

    /**
     * @brief Computes the first variation of the tangent derivative.
     */
    void ComputeTangentDerivativeVariation(
        const Vector& rShapeFunctionDerivatives,
        const Vector& rShapeFunctionSecondDerivatives,
        const array_1d<double, 3>& rFirstBaseVector,
        const array_1d<double, 3>& rSecondBaseVector,
        Vector& rTangentDerivativeVariation) const;

    /**
     * @brief Computes the second variation of the tangent derivative.
     */
    void ComputeTangentDerivativeSecondVariation(
        const Vector& rShapeFunctionDerivatives,
        const Vector& rShapeFunctionSecondDerivatives,
        const array_1d<double, 3>& rFirstBaseVector,
        const array_1d<double, 3>& rSecondBaseVector,
        Matrix& rTangentDerivativeSecondVariation) const;

    /**
     * @brief Computes the first variation of the second tangent derivative.
     */
    void ComputeTangentSecondDerivativeVariation(
        const Vector& rShapeFunctionDerivatives,
        const Vector& rShapeFunctionSecondDerivatives,
        const Vector& rShapeFunctionThirdDerivatives,
        const array_1d<double, 3>& rFirstBaseVector,
        const array_1d<double, 3>& rSecondBaseVector,
        const array_1d<double, 3>& rThirdBaseVector,
        Vector& rTangentSecondDerivativeVariation) const;

    /**
     * @brief Computes the Rodrigues rotation matrix.
     */
    static void ComputeRodriguesMatrix(
        const array_1d<double, 3>& rAxis,
        const double Phi,
        Matrix3d& rRodriguesMatrix);

    /**
     * @brief Computes the first derivative of the Rodrigues rotation matrix.
     */
    static void ComputeRodriguesMatrixDerivative(
        const array_1d<double, 3>& rAxis,
        const array_1d<double, 3>& rAxisDerivative,
        const double Phi,
        const double PhiDerivative,
        Matrix3d& rRodriguesMatrixDerivative);

    /**
     * @brief Computes the second derivative of the Rodrigues rotation matrix.
     */
    static void ComputeRodriguesMatrixSecondDerivative(
        const array_1d<double, 3>& rAxis,
        const array_1d<double, 3>& rAxisDerivative,
        const array_1d<double, 3>& rAxisSecondDerivative,
        const double Phi,
        const double PhiDerivative,
        const double PhiSecondDerivative,
        Matrix3d& rRodriguesMatrixSecondDerivative);

    static void ComputeRodriguesMatrixVariation(
        const array_1d<double, 3>& rAxis,
        const Vector& rAxisVariation,
        const Vector& rShapeFunctions,
        const double Phi,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rRodriguesMatrixVariation);

    static void ComputeRodriguesMatrixSecondVariation(
        const array_1d<double, 3>& rAxis,
        const Vector& rAxisVariation,
        const Matrix& rAxisSecondVariation,
        const Vector& rShapeFunctions,
        const double Phi,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rRodriguesMatrixSecondVariation);

    static void ComputeRodriguesMatrixDerivativeVariation(
        const array_1d<double, 3>& rAxis,
        const Vector& rAxisVariation,
        const array_1d<double, 3>& rAxisDerivative,
        const Vector& rAxisDerivativeVariation,
        const Vector& rShapeFunctions,
        const Vector& rShapeFunctionDerivatives,
        const double Phi,
        const double PhiDerivative,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rRodriguesMatrixDerivativeVariation);

    static void ComputeRodriguesMatrixDerivativeSecondVariation(
        const array_1d<double, 3>& rAxis,
        const Vector& rAxisVariation,
        const array_1d<double, 3>& rAxisDerivative,
        const Vector& rAxisDerivativeVariation,
        const Matrix& rAxisSecondVariation,
        const Matrix& rAxisDerivativeSecondVariation,
        const Vector& rShapeFunctions,
        const Vector& rShapeFunctionDerivatives,
        const double Phi,
        const double PhiDerivative,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rRodriguesMatrixDerivativeSecondVariation);

    static void ComputeRodriguesMatrixSecondDerivativeVariation(
        const array_1d<double, 3>& rAxis,
        const Vector& rAxisVariation,
        const array_1d<double, 3>& rAxisDerivative,
        const Vector& rAxisDerivativeVariation,
        const array_1d<double, 3>& rAxisSecondDerivative,
        const Vector& rAxisSecondDerivativeVariation,
        const Vector& rShapeFunctions,
        const Vector& rShapeFunctionDerivatives,
        const Vector& rShapeFunctionSecondDerivatives,
        const double Phi,
        const double PhiDerivative,
        const double PhiSecondDerivative,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rRodriguesMatrixSecondDerivativeVariation);

    static void ComputeRodriguesMatrixVariations(
        const array_1d<double, 3>& rAxis,
        const array_1d<double, 3>& rAxisVariation,
        const array_1d<double, 3>& rAxisDerivative,
        const Vector& rAxisDerivativeVariation,
        const Matrix& rAxisSecondVariation,
        const Matrix& rAxisDerivativeSecondVariation,
        const Vector& rShapeFunctions,
        const Vector& rShapeFunctionDerivatives,
        const double Phi,
        const double PhiDerivative,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rRodriguesMatrixVariation,
        Matrix& rRodriguesMatrixDerivativeVariation,
        Matrix& rRodriguesMatrixSecondVariation,
        Matrix& rRodriguesMatrixDerivativeSecondVariation);

    /**
     * @brief Computes the rotation mapping between reference and current tangents.
     */
    static void ComputeLambdaMatrix(
        const array_1d<double, 3>& rReferenceTangent,
        const array_1d<double, 3>& rCurrentTangent,
        Matrix3d& rLambdaMatrix);

    static void ComputeLambdaMatrixDerivative(
        const array_1d<double, 3>& rReferenceTangent,
        const array_1d<double, 3>& rCurrentTangent,
        const array_1d<double, 3>& rReferenceTangentDerivative,
        const array_1d<double, 3>& rCurrentTangentDerivative,
        Matrix3d& rLambdaMatrixDerivative);

    static void ComputeLambdaMatrixSecondDerivative(
        const array_1d<double, 3>& rReferenceTangent,
        const array_1d<double, 3>& rCurrentTangent,
        const array_1d<double, 3>& rReferenceTangentDerivative,
        const array_1d<double, 3>& rCurrentTangentDerivative,
        const array_1d<double, 3>& rReferenceTangentSecondDerivative,
        const array_1d<double, 3>& rCurrentTangentSecondDerivative,
        Matrix3d& rLambdaMatrixSecondDerivative);

    static void ComputeLambdaMatrixVariation(
        const array_1d<double, 3>& rReferenceTangent,
        const array_1d<double, 3>& rCurrentTangent,
        const Vector& rCurrentTangentVariation,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rLambdaMatrixVariation);

    static void ComputeLambdaMatrixSecondVariation(
        const array_1d<double, 3>& rReferenceTangent,
        const array_1d<double, 3>& rCurrentTangent,
        const Vector& rCurrentTangentVariation,
        const Matrix& rCurrentTangentSecondVariation,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rLambdaMatrixSecondVariation);

    static void ComputeLambdaMatrixDerivativeVariation(
        const array_1d<double, 3>& rReferenceTangent,
        const array_1d<double, 3>& rCurrentTangent,
        const array_1d<double, 3>& rReferenceTangentDerivative,
        const Vector& rCurrentTangentVariation,
        const array_1d<double, 3>& rCurrentTangentDerivative,
        const Vector& rCurrentTangentDerivativeVariation,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rLambdaMatrixDerivativeVariation);

    static void ComputeLambdaMatrixSecondDerivativeVariation(
        const array_1d<double, 3>& rReferenceTangent,
        const array_1d<double, 3>& rCurrentTangent,
        const array_1d<double, 3>& rReferenceTangentDerivative,
        const array_1d<double, 3>& rReferenceTangentSecondDerivative,
        const Vector& rCurrentTangentVariation,
        const array_1d<double, 3>& rCurrentTangentDerivative,
        const array_1d<double, 3>& rCurrentTangentSecondDerivative,
        const Vector& rCurrentTangentDerivativeVariation,
        const Vector& rCurrentTangentSecondDerivativeVariation,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rLambdaMatrixSecondDerivativeVariation);

    static void ComputeLambdaMatrixDerivativeSecondVariation(
        const array_1d<double, 3>& rReferenceTangent,
        const array_1d<double, 3>& rCurrentTangent,
        const array_1d<double, 3>& rReferenceTangentDerivative,
        const Vector& rCurrentTangentVariation,
        const array_1d<double, 3>& rCurrentTangentDerivative,
        const Vector& rCurrentTangentDerivativeVariation,
        const Matrix& rCurrentTangentSecondVariation,
        const Matrix& rCurrentTangentDerivativeSecondVariation,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rLambdaMatrixDerivativeSecondVariation);

    static void ComputeLambdaMatrixVariations(
        const array_1d<double, 3>& rReferenceTangent,
        const array_1d<double, 3>& rCurrentTangent,
        const array_1d<double, 3>& rReferenceTangentDerivative,
        const Vector& rCurrentTangentVariation,
        const array_1d<double, 3>& rCurrentTangentDerivative,
        const Vector& rCurrentTangentDerivativeVariation,
        const Matrix& rCurrentTangentSecondVariation,
        const Matrix& rCurrentTangentDerivativeSecondVariation,
        const SizeType NumberOfDofs,
        const SizeType DofsPerNode,
        Matrix& rLambdaMatrixVariation,
        Matrix& rLambdaMatrixDerivativeVariation,
        Matrix& rLambdaMatrixSecondVariation,
        Matrix& rLambdaMatrixDerivativeSecondVariation);

    ///@}
    ///@name Member Variables
    ///@{

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    SizeType mNumberOfDofs = 0;
    SizeType mDofsPerNode = 4;

    ///@}

private:

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(
        Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(
            rSerializer,
            Element);
    }

    void load(
        Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(
            rSerializer,
            Element);
    }

    ///@}

}; // Class BeamBaseElement3D

} // namespace Kratos