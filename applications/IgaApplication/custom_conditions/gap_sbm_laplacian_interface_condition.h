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
#include "includes/define.h"
#include "includes/condition.h"

// Project includes
#include "includes/variables.h"
#include "iga_application_variables.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) GapSbmLaplacianInterfaceCondition : public Condition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GapSbmLaplacianInterfaceCondition);

    GapSbmLaplacianInterfaceCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry) {}
    GapSbmLaplacianInterfaceCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties) {}
    GapSbmLaplacianInterfaceCondition() = default;
    ~GapSbmLaplacianInterfaceCondition() override = default;

    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GapSbmLaplacianInterfaceCondition>(NewId, pGeom, pProperties);
    }
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GapSbmLaplacianInterfaceCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;
    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"GapSbmLaplacianInterfaceCondition\" #" << Id();
        return buffer.str();
    }

protected:
    void InitializeMemberVariables();
    void InitializeSbmMemberVariables();

    void ComputeTaylorExpansionContribution(const GeometryType& rGeometry, const Vector& rDistanceVector, Vector& H_sum_vec);
    void ComputeGradientTaylorExpansionContribution(const GeometryType& rGeometry, const Vector& rDistanceVector, Matrix& grad_H_sum);
    double ComputeTaylorTerm(double derivative, double dx, IndexType n_k, double dy, IndexType k);
    double ComputeTaylorTerm3D(double derivative, double dx, IndexType k_x, double dy, IndexType k_y, double dz, IndexType k_z);

    const GeometryType& GetGeometryPlus() const { return *this->GetValue(NEIGHBOUR_GEOMETRIES)[0]; }
    const GeometryType& GetGeometryMinus() const { return *this->GetValue(NEIGHBOUR_GEOMETRIES)[1]; }

    // state
    unsigned int mDim = 0;
    IndexType mBasisFunctionsOrder = 0;
    array_1d<double,3> mNormalPhysicalSpace;
    Vector mDistanceVectorPlus;
    Vector mDistanceVectorMinus;
    double mPenalty = 0.0; // scaled by p^2 / h
    double mNitschePenalty = 1.0; // 1.0 penalty, -1.0 penalty-free

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override { KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition); }
    void load(Serializer& rSerializer) override { KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition); }
};

} // namespace Kratos

