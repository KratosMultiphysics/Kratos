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
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"

// Application includes
#include "iga_application_variables.h"
#include "custom_elements/beam_base_element_3D.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class NonLinearBernoulliBeamElement3D
 * @brief Three-dimensional Bernoulli beam element.
 */
class KRATOS_API(IGA_APPLICATION) NonLinearBernoulliBeamElement3D
    : public BeamBaseElement3D
{
protected:

    using Vector3d = array_1d<double, 3>;
    using Matrix3d = BoundedMatrix<double, 3, 3>;
    Matrix mSMatrixRodriguesVariation;
    Matrix mSMatrixLambdaVariation;
    Matrix mSMatrixLambdaVariationRodriguesLambda;
    Matrix mSMatrixRodriguesLambdaVariationRodriguesLambda;
    Matrix mSMatrixRodriguesVariationLambdaRodriguesLambda;
    Matrix mSMatrixRodriguesLambdaRodriguesLambdaVariation;

    Matrix mSMatrixRodriguesDerivativeVariation;
    Matrix mSMatrixLambdaDerivativeVariation;
    Matrix mSMatrixLambdaDerivativeVariationRodriguesLambda;
    Matrix mSMatrixLambdaVariationRodriguesDerivativeLambda;
    Matrix mSMatrixLambdaVariationRodriguesLambdaDerivative;

    Matrix mSMatrixRodriguesDerivativeLambdaVariationRodriguesLambda;
    Matrix mSMatrixRodriguesLambdaDerivativeVariationRodriguesLambda;
    Matrix mSMatrixRodriguesLambdaVariationRodriguesDerivativeLambda;
    Matrix mSMatrixRodriguesLambdaVariationRodriguesLambdaDerivative;

    Matrix mSMatrixRodriguesDerivativeVariationLambdaRodriguesLambda;
    Matrix mSMatrixRodriguesVariationLambdaDerivativeRodriguesLambda;
    Matrix mSMatrixRodriguesVariationLambdaRodriguesDerivativeLambda;
    Matrix mSMatrixRodriguesVariationLambdaRodriguesLambdaDerivative;
    Matrix mSMatrixRodriguesLambdaRodriguesLambdaDerivativeVariation;

    Matrix mSMatrixRodriguesSecondVariation;
    Matrix mSMatrixLambdaSecondVariation;
    Matrix mSMatrixLambdaSecondVariationRodriguesLambda;
    Matrix mSMatrixRodriguesDerivativeSecondVariation;
    Matrix mSMatrixLambdaDerivativeSecondVariation;
    Matrix mSMatrixLambdaDerivativeSecondVariationRodriguesLambda;
    Matrix mSMatrixLambdaSecondVariationRodriguesDerivativeLambda;
    Matrix mSMatrixLambdaSecondVariationRodriguesLambdaDerivative;

    struct KinematicVariables
    {
        // Reference configuration
        Vector3d R1;
        Vector3d R2;
        Vector3d R3;

        double A;
        double B;

        // Current configuration
        Vector3d r1;
        Vector3d r2;
        Vector3d r3;

        double a;
        double b;

        // Cross-section directors
        Vector3d N0;
        Vector3d V0;

        Vector3d n;
        Vector3d v;

        // Curvature components
        double B_n;
        double B_v;

        double b_n;
        double b_v;

        double C_12;
        double C_13;

        double c_12;
        double c_13;

        // Rotations
        double Phi;
        double Phi_der;

        double phi;
        double phi_der;

        explicit KinematicVariables(
            const SizeType Dimension = 3)
        {
            R1 = ZeroVector(Dimension);
            R2 = ZeroVector(Dimension);
            R3 = ZeroVector(Dimension);

            r1 = ZeroVector(Dimension);
            r2 = ZeroVector(Dimension);
            r3 = ZeroVector(Dimension);

            N0 = ZeroVector(Dimension);
            V0 = ZeroVector(Dimension);

            n = ZeroVector(Dimension);
            v = ZeroVector(Dimension);

            A = 0.0;
            B = 0.0;

            a = 0.0;
            b = 0.0;

            B_n = 0.0;
            B_v = 0.0;

            b_n = 0.0;
            b_v = 0.0;

            C_12 = 0.0;
            C_13 = 0.0;

            c_12 = 0.0;
            c_13 = 0.0;

            Phi = 0.0;
            Phi_der = 0.0;

            phi = 0.0;
            phi_der = 0.0;
        }

    };

public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(
        NonLinearBernoulliBeamElement3D);

    using BaseType = BeamBaseElement3D;
    using SizeType = typename BaseType::SizeType;
    using IndexType = typename BaseType::IndexType;
    using GeometryType = typename BaseType::GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    NonLinearBernoulliBeamElement3D(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : BeamBaseElement3D(
            NewId,
            pGeometry)
    {
    }

    NonLinearBernoulliBeamElement3D(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : BeamBaseElement3D(
            NewId,
            pGeometry,
            pProperties)
    {
    }

    NonLinearBernoulliBeamElement3D()
        : BeamBaseElement3D()
    {
    }

    ~NonLinearBernoulliBeamElement3D() override = default;

    ///@}
    ///@name Element Creation
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<NonLinearBernoulliBeamElement3D>(
            NewId,
            pGeometry,
            pProperties);
    }

    Element::Pointer Create(
        IndexType NewId,
        const NodesArrayType& rThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<NonLinearBernoulliBeamElement3D>(
            NewId,
            GetGeometry().Create(rThisNodes),
            pProperties);
    }

    ///@}
    ///@name Element Calculations
    ///@{

    void Initialize(
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateKinematics(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables);

    void CalculateConstitutiveVariables(
        const IndexType IntegrationPointIndex,
        const KinematicVariables& rKinematicVariables,
        ConstitutiveVariables& rConstitutiveVariables,
        ConstitutiveLaw::Parameters& rConstitutiveLawParameters,
        const ConstitutiveLaw::StressMeasure StressMeasure);
    
    void ComputeBMatrices(
        const IndexType IntegrationPointIndex,
        const KinematicVariables& rKinematicVariables,
        Matrix& rBAxial,
        Matrix& rBBending1,
        Matrix& rBBending2,
        Matrix& rBTorsion1,
        Matrix& rBTorsion2);

    void ComputeGMatrices(
        const IndexType IntegrationPointIndex,
        const KinematicVariables& rKinematicVariables,
        Matrix& rGAxial,
        Matrix& rGBending1,
        Matrix& rGBending2,
        Matrix& rGTorsion1,
        Matrix& rGTorsion2);
    
        void ComputeGAxial(
        const IndexType IntegrationPointIndex,
        Matrix& rGAxial,
        KinematicVariables& rKinematicVariables) const;

    void ComputeGBending(
        const IndexType IntegrationPointIndex,
        Matrix& rGBending1,
        Matrix& rGBending2,
        KinematicVariables& rKinematicVariables) const;

    void ComputeGTorsion(
        const IndexType IntegrationPointIndex,
        Matrix& rGTorsion1,
        Matrix& rGTorsion2,
        KinematicVariables& rKinematicVariables) const;

    void ComputeBAxial(
        const IndexType IntegrationPointIndex,
        Matrix& rBAxial,
        KinematicVariables& rKinematicVariables) const;

    void ComputeBBending(
        const IndexType IntegrationPointIndex,
        Matrix& rBBending1,
        Matrix& rBBending2,
        KinematicVariables& rKinematicVariables) const;

    void ComputeBTorsion(
        const IndexType IntegrationPointIndex,
        Matrix& rBTorsion1,
        Matrix& rBTorsion2,
        KinematicVariables& rKinematicVariables) const;

    ///@}
    ///@name Input and Output
    ///@{

    int Check(
        const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}

private:

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType number_of_nodes =
            GetGeometry().size();

        const SizeType number_of_dofs =
            number_of_nodes * mDofsPerNode;

        if (rRightHandSideVector.size() != number_of_dofs) {
            rRightHandSideVector.resize(
                number_of_dofs,
                false);
        }

        noalias(rRightHandSideVector) =
            ZeroVector(number_of_dofs);

        MatrixType left_hand_side_matrix;

        CalculateAll(
            left_hand_side_matrix,
            rRightHandSideVector,
            rCurrentProcessInfo,
            false,
            true);
    }

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType number_of_nodes =
            GetGeometry().size();

        const SizeType number_of_dofs =
            number_of_nodes * mDofsPerNode;

        VectorType right_hand_side_vector;

        if (
            rLeftHandSideMatrix.size1() != number_of_dofs ||
            rLeftHandSideMatrix.size2() != number_of_dofs)
        {
            rLeftHandSideMatrix.resize(
                number_of_dofs,
                number_of_dofs,
                false);
        }

        noalias(rLeftHandSideMatrix) =
            ZeroMatrix(
                number_of_dofs,
                number_of_dofs);

        CalculateAll(
            rLeftHandSideMatrix,
            right_hand_side_vector,
            rCurrentProcessInfo,
            true,
            false);
    }

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType number_of_nodes =
            GetGeometry().size();

        const SizeType number_of_dofs =
            number_of_nodes * mDofsPerNode;

        if (rRightHandSideVector.size() != number_of_dofs) {
            rRightHandSideVector.resize(
                number_of_dofs,
                false);
        }

        noalias(rRightHandSideVector) =
            ZeroVector(number_of_dofs);

        if (
            rLeftHandSideMatrix.size1() != number_of_dofs ||
            rLeftHandSideMatrix.size2() != number_of_dofs)
        {
            rLeftHandSideMatrix.resize(
                number_of_dofs,
                number_of_dofs,
                false);
        }

        noalias(rLeftHandSideMatrix) =
            ZeroMatrix(
                number_of_dofs,
                number_of_dofs);

        CalculateAll(
            rLeftHandSideMatrix,
            rRightHandSideVector,
            rCurrentProcessInfo,
            true,
            true);
    }

    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag) override;

    void ComputeReferenceTwistAngleAndDerivative(
        const KinematicVariables& rKinematicVariables,
        double& rPhi,
        double& rPhiDerivative);

    double CalculateDeltaPhi(
        const KinematicVariables& rKinematicVariables,
        const Vector3d& rNormal);

    void ComputeReferenceCrossSectionGeometry(
        KinematicVariables& rKinematicVariables);

    void ComputeCurrentCrossSectionGeometry(
        KinematicVariables& rKinematicVariables);

    friend class Serializer;

    void save(
        Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(
            rSerializer,
            BeamBaseElement3D);
    }

    void load(
        Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(
            rSerializer,
            BeamBaseElement3D);
    }
};

///@}

} // Namespace Kratos