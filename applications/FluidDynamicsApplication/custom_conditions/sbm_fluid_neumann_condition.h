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

// Project includes
#include "fluid_dynamics_application_variables.h"
#include "includes/condition.h"

namespace Kratos {

/**
 * @class SbmFluidNeumannCondition2D4N
 * @brief Applies a prescribed traction on a reconstructed Gap-SBM boundary.
 * @details The condition geometry is the owner Quad4, so the extended Q1
 *          basis is evaluated directly at the two reconstructed-boundary
 *          Gauss points stored in SURROGATE_BOUNDARY_PROJECTION. The 2x7
 *          matrix rows contain [x,y,z,nx,ny,nz,weight]. A pointwise prescribed
 *          total Cauchy stress can be stored in CAUCHY_STRESS_TENSOR as a 2x4
 *          matrix [sigma_xx,sigma_xy,sigma_yx,sigma_yy]; the condition then
 *          applies t=sigma*n. Alternatively, FACE_LOAD supplies one constant
 *          traction vector for the complete face.
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) SbmFluidNeumannCondition2D4N
    : public Condition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SbmFluidNeumannCondition2D4N);

    using IndexType = Condition::IndexType;
    using SizeType = Condition::SizeType;
    using GeometryType = Condition::GeometryType;
    using NodesArrayType = Condition::NodesArrayType;
    using PropertiesType = Condition::PropertiesType;
    using MatrixType = Condition::MatrixType;
    using VectorType = Condition::VectorType;
    using EquationIdVectorType = Condition::EquationIdVectorType;
    using DofsVectorType = Condition::DofsVectorType;

    SbmFluidNeumannCondition2D4N(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    SbmFluidNeumannCondition2D4N(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    ~SbmFluidNeumannCondition2D4N() override = default;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) const override;

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesType::Pointer pProperties) const override;

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
        return "SbmFluidNeumannCondition2D4N";
    }

protected:
    SbmFluidNeumannCondition2D4N() = default;

private:
    friend class Serializer;

    void CalculateNeumannLoad(VectorType& rRightHandSideVector) const;

    array_1d<double, 3> GetPrescribedTraction(
        SizeType IntegrationPointIndex,
        const array_1d<double, 2>& rNormal) const;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }
};

} // namespace Kratos
