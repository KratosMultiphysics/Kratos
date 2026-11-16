// BSD License, Kratos default license: kratos/license.txt
#pragma once

#include "includes/condition.h"

namespace Kratos {

/// Symmetric Nitsche Dirichlet condition on the physical Gap-SBM boundary.
/// NEIGHBOUR_ELEMENTS stores one Q1 parent, and condition geometry contains
/// its nodes. GAP_INTERFACE_COORDINATES stores the oriented physical
/// segment. GAP_DIRICHLET_VALUES contains three values, ordered as the
/// line Gauss-3 points returned by CalculateOnIntegrationPoints.
/// The implemented flux and penalty assume unit diffusivity and do not
/// include advective fluxes.
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) LaplacianGapDirichletCondition2D : public Condition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LaplacianGapDirichletCondition2D);
    using Condition::Condition;
    Condition::Pointer Create(IndexType Id, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const override;
    Condition::Pointer Create(IndexType Id, const NodesArrayType& rNodes, PropertiesType::Pointer pProperties) const override;
    void CalculateLocalSystem(MatrixType& rLHS, VectorType& rRHS, const ProcessInfo& rInfo) override;
    using Condition::CalculateOnIntegrationPoints;
    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rValues, const ProcessInfo& rInfo) override;
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rInfo) const override;
    void GetDofList(DofsVectorType& rDofs, const ProcessInfo& rInfo) const override;
    int Check(const ProcessInfo& rInfo) const override;
protected:
    virtual bool IsConnector() const { return false; }
};

/// Symmetric interior-penalty coupling between two extended parent fields.
/// NEIGHBOUR_ELEMENTS stores the two parents and the condition geometry is
/// their node union. SURROGATE_BOUNDARY_FACE_COORDINATES stores the oriented
/// connector; its right normal points from the first parent to the second.
/// This condition also assumes unit diffusivity and has no advective flux.
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) LaplacianGapConnectorCondition2D : public LaplacianGapDirichletCondition2D
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LaplacianGapConnectorCondition2D);
    using LaplacianGapDirichletCondition2D::LaplacianGapDirichletCondition2D;
    Condition::Pointer Create(IndexType Id, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const override;
    Condition::Pointer Create(IndexType Id, const NodesArrayType& rNodes, PropertiesType::Pointer pProperties) const override;
protected:
    bool IsConnector() const override { return true; }
};
}
