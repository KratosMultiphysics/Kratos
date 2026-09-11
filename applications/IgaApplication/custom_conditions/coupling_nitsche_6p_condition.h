//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//

#pragma once

// System includes
#include "includes/define.h"
#include "includes/condition.h"

// External includes

// Project includes
#include "iga_application_variables.h"
#include "custom_utilities/iga_flags.h"

#include "geometries/coupling_geometry.h"

namespace Kratos
{

/// Nitsche coupling condition for the degenerated-solid Reissner-Mindlin
/// shell_6p_element.
/**
*   Enforces displacement AND rotation continuity between two shell_6p
*   patches through a single, physically consistent Nitsche term
*/
class CouplingNitsche6pCondition
    : public Condition
{
protected:

    ///@name Internal structs
    ///@{

    /// Reference-configuration midsurface kinematics at one boundary
    /// integration point, for one patch (master or slave).
    struct KinematicVariables
    {
        array_1d<double, 3> BaseVector1;
        array_1d<double, 3> BaseVector2;
        array_1d<double, 3> NormalVector;      // a3, normalized
        array_1d<double, 3> NormalVectorTilde;  // a1 x a2, not normalized
        double DifferentialArea;

        explicit KinematicVariables(IndexType Dimension)
        {
            noalias(BaseVector1) = ZeroVector(Dimension);
            noalias(BaseVector2) = ZeroVector(Dimension);
            noalias(NormalVector) = ZeroVector(Dimension);
            noalias(NormalVectorTilde) = ZeroVector(Dimension);
            DifferentialArea = 1.0;
        }
    };

    ///@}

public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of CouplingNitsche6pCondition
    KRATOS_CLASS_POINTER_DEFINITION(CouplingNitsche6pCondition);

    /// Index type
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry
    CouplingNitsche6pCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {};

    /// Constructor with Id, geometry and property
    CouplingNitsche6pCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {};

    /// Default constructor
    CouplingNitsche6pCondition()
        : Condition()
    {};

    /// Destructor.
    virtual ~CouplingNitsche6pCondition() = default;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Create with Id, pointer to geometry and pointer to property
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<CouplingNitsche6pCondition>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive< CouplingNitsche6pCondition >(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    ///@}
    ///@name Operations
    ///@{

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        MatrixType left_hand_side_matrix = Matrix(0, 0);

        CalculateAll(left_hand_side_matrix, rRightHandSideVector,
            rCurrentProcessInfo, false, true);
    }

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        VectorType right_hand_side_vector = Vector(0);

        CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
            rCurrentProcessInfo, true, false);
    }

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
            rCurrentProcessInfo, true, true);
    }

    /**
    * @brief Sets on rResult the ID's of the element degrees of freedom
    */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    /**
    * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
    */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    /// Calculates left (K) and right (u) hand sides, according to the flags
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    );

    ///@}
    ///@name Check
    ///@{

    /// Performs check if NITSCHE_STABILIZATION_FACTOR is provided.
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"CouplingNitsche6pCondition\" #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"CouplingNitsche6pCondition\" #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}

private:

    ///@name Private enums
    ///@{

    enum class PatchType { Master, Slave };

    ///@}
    ///@name Private operations
    ///@{

    /// Reference-config midsurface kinematics (a1, a2, a3, a3_tilde, dA) at
    /// one boundary integration point of one patch.
    void CalculateKinematics(
        IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        const PatchType& rPatch) const;

    /// d(a3)/d(xi), d(a3)/d(eta) 
    void CalculateNormalVectorDerivatives(
        IndexType IntegrationPointIndex,
        const KinematicVariables& rKinematicVariables,
        Matrix& rNormalVectorDerivatives,
        const PatchType& rPatch) const;

    /// Voigt transformation matrix from the local (a1,a2,a3) orthonormal
    /// frame to global Cartesian
    void CalculateTransformationFromLocalToGlobalCartesian(
        const KinematicVariables& rKinematicVariables,
        Matrix& rTransformationMatrix) const;

    /// Local orthotropic-in-plane / transverse-shear linear elastic
    /// constitutive matrix (6x6, Voigt [xx,yy,zz,xy,yz,xz])
    void CalculateLocalConstitutiveMatrix(
        const Properties& rProperties,
        Matrix& rConstitutiveMatrixLocal) const;

    /// Physical lateral-surface conormal (unit vector, in the tangent plane)
    /// via Nanson's formula: n_phys = det(J) * J^-T * nu_param
    void CalculateLateralConormal(
        IndexType IntegrationPointIndex,
        const Matrix& rJacobianInv,
        double JacobianThicknessDet,
        const PatchType& rPatch,
        array_1d<double, 3>& rUnitConormal,
        double& rAreaScale) const;

    /// Builds the 3x(6*n_nodes) global-Cartesian strain B-operator at one
    /// through-thickness point (zeta)
    void CalculateBAndDisplacementOperator(
        IndexType IntegrationPointIndex,
        double zeta,
        double Thickness,
        const Matrix& rJacobianInv,
        const Matrix& rNormalVectorDerivatives,
        const KinematicVariables& rKinematicVariables,
        const PatchType& rPatch,
        Matrix& rBOperator,
        Matrix& rDisplacementOperator) const;

    /// Assembles the 3x3 thickness Jacobian at one (surface IP, zeta) point
    /// and its inverse/determinant
    void CalculateThicknessJacobian(
        double zeta,
        double Thickness,
        const KinematicVariables& rKinematicVariables,
        const Matrix& rNormalVectorDerivatives,
        Matrix& rJacobianInv,
        double& rJacobianDet) const;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }

    ///@}

}; // Class CouplingNitsche6pCondition

}  // namespace Kratos.
