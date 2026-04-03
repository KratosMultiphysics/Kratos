// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    GitHub Copilot
//

#if !defined(KRATOS_CONSISTENT_FLUX_BOUNDARY_CONDITION_H_INCLUDED)
#define KRATOS_CONSISTENT_FLUX_BOUNDARY_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/global_pointers_vector.h"
#include "geometries/geometry.h"
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"

namespace Kratos
{

class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) ConsistentFluxBoundaryCondition : public Condition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ConsistentFluxBoundaryCondition);

    typedef Condition::MatrixType MatrixType;
    typedef Condition::VectorType VectorType;

    ConsistentFluxBoundaryCondition(IndexType NewId, Geometry<Node>::Pointer pGeometry);

    ConsistentFluxBoundaryCondition(
        IndexType NewId,
        Geometry<Node>::Pointer pGeometry,
        Properties::Pointer pProperties);

    ~ConsistentFluxBoundaryCondition() override;

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

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

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& rConditionalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    std::string Info() const override;

    void PrintInfo(std::ostream& rOStream) const override;

    void PrintData(std::ostream& rOStream) const override;

protected:
    ConsistentFluxBoundaryCondition();

private:
    struct IntegrationPointData
    {
        double Weight = 0.0;
        double Conductivity = 0.0;
        Vector ParentN;
        Matrix ParentDNDX;
        array_1d<double, 3> UnitNormal = array_1d<double, 3>{0.0, 0.0, 0.0};
    };

    Element& GetParentElement();

    const Element& GetParentElement() const;

    void CalculateConditionSystem(
        MatrixType* pLeftHandSideMatrix,
        VectorType* pRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) const;

    void ComputeIntegrationPointData(
        const GeometryType& rConditionGeometry,
        const GeometryType& rParentGeometry,
        const IndexType IntegrationPointIndex,
        const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
        const Vector& rDetJ,
        const Variable<double>& rDiffusionVariable,
        IntegrationPointData& rData) const;

    void CalculateParentGeometryValues(
        const GeometryType& rConditionGeometry,
        const GeometryType& rParentGeometry,
        const GeometryType::CoordinatesArrayType& rConditionLocalCoordinates,
        Vector& rParentN,
        Matrix& rParentDNDX) const;

    Vector GetParentUnknownValues(const Variable<double>& rUnknownVariable) const;

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ConsistentFluxBoundaryCondition& operator=(ConsistentFluxBoundaryCondition const& rOther);

    ConsistentFluxBoundaryCondition(ConsistentFluxBoundaryCondition const& rOther);
};

} // namespace Kratos

#endif // KRATOS_CONSISTENT_FLUX_BOUNDARY_CONDITION_H_INCLUDED