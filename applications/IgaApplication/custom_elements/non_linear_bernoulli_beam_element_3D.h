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
 * @brief Geometrically nonlinear isogeometric spatial Bernoulli beam with torsion.
 * @details Total Lagrangian formulation derived from the 3D continuum with Green-Lagrange strains and
 * second Piola-Kirchhoff stresses. The cross section remains orthogonal to the center line (Bernoulli)
 * and warping is neglected. Each control point has four dofs: DISPLACEMENT_X/Y/Z and the rotation
 * about the center line (ROTATION_X).
 *
 * The cross-section frame is A_alpha = R_T(Psi) Lambda(T0, T) A0_alpha in the reference configuration and
 * a_alpha = R_t(psi) Lambda(T, t) A_alpha in the current one, where Lambda maps one tangent onto the other
 * and R rotates about the tangent. The generalized strains are the axial strain, the two bending curvatures
 * and the two torsional curvatures.
 *
 * Reference:
 * A.M. Bauer, M. Breitenberger, B. Philipp, R. Wuechner, K.-U. Bletzinger,
 * "Nonlinear isogeometric spatial Bernoulli beam",
 * Computer Methods in Applied Mechanics and Engineering 303 (2016) 101-127.
 * https://doi.org/10.1016/j.cma.2015.12.027
 * Equation numbers in the comments of this element refer to this paper.
 */
class KRATOS_API(IGA_APPLICATION) NonLinearBernoulliBeamElement3D
    : public BeamBaseElement3D
{
protected:

    using Vector3d = array_1d<double, 3>;
    using Matrix3d = BoundedMatrix<double, 3, 3>;

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

    /// Velocities [VELOCITY_X, VELOCITY_Y, VELOCITY_Z, ANGULAR_VELOCITY_X] per control point.
    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0) const override;

    /// Accelerations [ACCELERATION_X, ACCELERATION_Y, ACCELERATION_Z, ANGULAR_ACCELERATION_X] per control point.
    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0) const override;

    /**
     * @brief Consistent mass matrix in the reference configuration.
     * @details M_ij = int N_i N_j rho A dL for the displacements and int N_i N_j rho (I22 + I33) dL
     * for the rotation about the center line (polar mass moment of inertia per unit length).
     */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
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

    /// Normalized initial tangent T_0, taken from the element data if set, otherwise from the properties.
    Vector3d GetInitialTangent() const;

    /// Initial normal N_0, taken from the element data if set, otherwise from the properties.
    Vector3d GetInitialNormal() const;

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
