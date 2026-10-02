//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

#pragma once

// System includes
#include <vector>

// Project includes
#include "fluid_dynamics_application_variables.h"
#include "includes/condition.h"
#include "includes/constitutive_law.h"
#include "includes/element.h"

namespace Kratos {

/**
 * @class SbmFluidGapInterfaceCondition2D
 * @brief Weakly couples the FEM extensions of two adjacent 2D gap patches.
 * @details The condition implements a nonsymmetric Nitsche interface form for
 *          the incompressible Stokes operator,
 *
 *          - [v] . {sigma(u,p)n} + {sigma(v,q)n} . [u]
 *          + gamma/h [v] . [u].
 *
 *          The condition geometry is the union of the nodes of both parent
 *          Quad4 elements. The two parents are supplied in NEIGHBOUR_ELEMENTS
 *          and the oriented connector (surrogate vertex to reconstructed skin)
 *          is supplied in SURROGATE_BOUNDARY_FACE_COORDINATES. The normal
 *          encoded by the connector orientation points from parent 0 to parent
 *          1. PENALTY_COEFFICIENT=0 selects the penalty-free form.
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION)
    SbmFluidGapInterfaceCondition2D : public Condition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(
        SbmFluidGapInterfaceCondition2D);

    using IndexType = Condition::IndexType;
    using SizeType = Condition::SizeType;
    using GeometryType = Condition::GeometryType;
    using NodesArrayType = Condition::NodesArrayType;
    using PropertiesType = Condition::PropertiesType;
    using MatrixType = Condition::MatrixType;
    using VectorType = Condition::VectorType;
    using EquationIdVectorType = Condition::EquationIdVectorType;
    using DofsVectorType = Condition::DofsVectorType;

    SbmFluidGapInterfaceCondition2D(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    SbmFluidGapInterfaceCondition2D(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    ~SbmFluidGapInterfaceCondition2D() override = default;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) const override;

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesType::Pointer pProperties) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalVelocityContribution(
        MatrixType& rDampingMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& rConditionDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    std::string Info() const override
    {
        return "SbmFluidGapInterfaceCondition2D";
    }

protected:
    SbmFluidGapInterfaceCondition2D() = default;

private:
    friend class Serializer;

    struct ConstitutiveVariables
    {
        ConstitutiveLaw::StrainVectorType StrainVector = ZeroVector(3);
        ConstitutiveLaw::StressVectorType StressVector = ZeroVector(3);
        ConstitutiveLaw::VoigtSizeMatrixType ConstitutiveMatrix =
            ZeroMatrix(3, 3);
    };

    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag);

    static void CalculateParentShapeFunctions(
        const GeometryType& rParentGeometry,
        const Point& rGlobalPoint,
        Vector& rN,
        Matrix& rDN_DX);

    static void CalculateB(Matrix& rB, const Matrix& rDN_DX);

    static void BuildStressFromVoigtColumn(
        Matrix& rStressTensor,
        const Matrix& rConstitutiveBMatrix,
        IndexType Column);

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLaws;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
        rSerializer.save("ConstitutiveLaws", mConstitutiveLaws);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
        rSerializer.load("ConstitutiveLaws", mConstitutiveLaws);
    }
};

} // namespace Kratos
