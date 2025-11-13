//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli

#pragma once

// System includes

// External includes

// Project includes
#include "includes/condition.h"
#include "includes/convection_diffusion_settings.h"
#include "iga_application_variables.h"

namespace Kratos
{

/**
 * @brief Shifted-boundary interface condition for Laplacian fields.
 *
 * The condition weakly enforces continuity of the unknown and diffusive flux
 * between two surrogate interfaces ("plus" and "minus") that approximate a
 * cut boundary. Taylor-expanded neighbour shape functions are employed to
 * evaluate the solution at the true boundary location.
 */
class KRATOS_API(IGA_APPLICATION) LaplacianCouplingCondition : public Condition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LaplacianCouplingCondition);

    using SizeType = std::size_t;
    using IndexType = std::size_t;

    LaplacianCouplingCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    LaplacianCouplingCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    LaplacianCouplingCondition() = default;
    ~LaplacianCouplingCondition() override = default;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<LaplacianCouplingCondition>(NewId, pGeom, pProperties);
    }

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<LaplacianCouplingCondition>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

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

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"LaplacianCouplingCondition\" #" << Id();
        return buffer.str();
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"LaplacianCouplingCondition\" #" << Id();
    }

    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    IntegrationMethod GetIntegrationMethod() const override
    {
        return GeometryData::IntegrationMethod::GI_GAUSS_1;
    }

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }

    void InitializeMemberVariables();

    void GetSolutionCoefficientVectorA(
        const Variable<double>& rUnknown,
        Vector& rValues) const;

    void GetSolutionCoefficientVectorB(
        const Variable<double>& rUnknown,
        Vector& rValues) const;

    const GeometryType& GetGeometryMirror() const
    {
        return *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    }
    
    array_1d<double, 3> mNormalParameterSpaceA = ZeroVector(3);
    array_1d<double, 3> mNormalParameterSpaceB = ZeroVector(3);

    array_1d<double, 3> mNormalPhysicalSpaceA = ZeroVector(3);
    array_1d<double, 3> mNormalPhysicalSpaceB = ZeroVector(3);

    void ComputeTaylorExpansionContribution(const GeometryType& rGeometry, const Vector& rDistanceVector, Vector& H_sum_vec);
    void ComputeGradientTaylorExpansionContribution(const GeometryType& rGeometry, const Vector& rDistanceVector, Matrix& grad_H_sum);
    double ComputeTaylorTerm(double derivative, double dx, IndexType n_k, double dy, IndexType k);
    double ComputeTaylorTerm3D(double derivative, double dx, IndexType k_x, double dy, IndexType k_y, double dz, IndexType k_z);

    unsigned int mDim = 0;
    IndexType mBasisFunctionsOrder = 0;
    bool mIsGapSbmCoupling = false;
    Vector mDistanceVectorB;
};

} // namespace Kratos
