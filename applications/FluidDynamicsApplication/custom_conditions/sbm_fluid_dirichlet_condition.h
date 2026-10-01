//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


#pragma once

// Project includes
#include "fluid_dynamics_application_variables.h"
#include "includes/condition.h"
#include "includes/constitutive_law.h"

namespace Kratos {

/**
 * @class SbmFluidDirichletCondition2D4N
 * @brief Face-based 2D SBM velocity Dirichlet condition.
 * @details One condition integrates an entire surrogate Line2 face while its
 *          Kratos geometry is the owner Quad4. This is required by standard
 *          schemes, which assume that geometry nodes and local-system blocks
 *          coincide. The two oriented face points are stored row-wise in
 *          SURROGATE_BOUNDARY_FACE_COORDINATES and true-boundary projections
 *          in SURROGATE_BOUNDARY_PROJECTION. Prescribed values at the face
 *          Gauss points are stored in SBM_BOUNDARY_VELOCITIES. The existing
 *          core PENALTY_COEFFICIENT controls the optional penalty contribution
 *          and defaults to zero (penalty-free).
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) SbmFluidDirichletCondition2D4N
    : public Condition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SbmFluidDirichletCondition2D4N);

    using IndexType = Condition::IndexType;
    using SizeType = Condition::SizeType;
    using GeometryType = Condition::GeometryType;
    using NodesArrayType = Condition::NodesArrayType;
    using PropertiesType = Condition::PropertiesType;
    using MatrixType = Condition::MatrixType;
    using VectorType = Condition::VectorType;
    using EquationIdVectorType = Condition::EquationIdVectorType;
    using DofsVectorType = Condition::DofsVectorType;

    SbmFluidDirichletCondition2D4N(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    SbmFluidDirichletCondition2D4N(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    ~SbmFluidDirichletCondition2D4N() override = default;

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
        return "SbmFluidDirichletCondition2D4N";
    }

protected:
    SbmFluidDirichletCondition2D4N() = default;

private:
    friend class Serializer;

    struct ConstitutiveVariables
    {
        ConstitutiveLaw::StrainVectorType StrainVector = ZeroVector(3);
        ConstitutiveLaw::StressVectorType StressVector = ZeroVector(3);
        ConstitutiveLaw::VoigtSizeMatrixType ConstitutiveMatrix = ZeroMatrix(3, 3);
    };

    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag);

    const GeometryType& GetParentGeometry() const;

    void CalculateParentShapeFunctions(
        const Point& rGlobalPoint,
        Vector& rN,
        Matrix& rDN_DX) const;

    array_1d<double, 3> GetPrescribedVelocity(SizeType IntegrationPointIndex) const;

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
