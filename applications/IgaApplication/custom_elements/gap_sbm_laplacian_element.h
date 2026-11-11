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
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"

// Application includes
#include "iga_application_variables.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) GapSbmLaplacianElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GapSbmLaplacianElement);

    GapSbmLaplacianElement(IndexType NewId, GeometryType::Pointer pGeometry);
    GapSbmLaplacianElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);
    GapSbmLaplacianElement() : Element() {}
    ~GapSbmLaplacianElement() override;

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;
    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    IntegrationMethod GetIntegrationMethod() const override;

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "GapSbmLaplacianElement #" << Id();
        return buffer.str();
    }

    // GAP-SBM helpers
    void InitializeMemberVariables();
    void InitializeSbmMemberVariables();

protected:
    // Taylor expansion helpers (copied pattern from GAP-SBM elements)
    void ComputeTaylorExpansionContribution(Vector& H_sum_vec);
    void ComputeGradientTaylorExpansionContribution(Matrix& grad_H_sum);
    double ComputeTaylorTerm(double derivative, double dx, IndexType n_k, double dy, IndexType k);
    double ComputeTaylorTerm3D(double derivative, double dx, IndexType k_x, double dy, IndexType k_y, double dz, IndexType k_z);

    const GeometryType& GetSurrogateGeometry() const
    {
        return *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    }

    std::size_t mDim = 0;
    Vector mDistanceVector; // between physical and surrogate centers
    IndexType mBasisFunctionsOrder = 0;

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }
};

} // namespace Kratos

