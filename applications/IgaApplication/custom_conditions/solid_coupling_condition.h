//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi

#pragma once

// System includes

// External includes

// Project includes
#include "includes/condition.h"
#include "iga_application_variables.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

/**
 * @brief Shifted-boundary interface condition for solid displacement fields.
 *
 * The condition weakly enforces continuity of the displacement and its normal
 * derivative between two surrogate interfaces ("plus" and "minus") that
 * approximate a cut boundary. Taylor-expanded neighbour shape functions are
 * employed to evaluate the solution at the true boundary location.
 */
class KRATOS_API(IGA_APPLICATION) SolidCouplingCondition : public Condition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SolidCouplingCondition);

    using SizeType = std::size_t;
    using IndexType = std::size_t;

    SolidCouplingCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    SolidCouplingCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    SolidCouplingCondition() = default;
    ~SolidCouplingCondition() override = default;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<SolidCouplingCondition>(NewId, pGeom, pProperties);
    }

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<SolidCouplingCondition>(
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
        buffer << "\"SolidCouplingCondition\" #" << Id();
        return buffer.str();
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"SolidCouplingCondition\" #" << Id();
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
    void InitializeMaterial();
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void GetSolutionCoefficientVectorA(Vector& rValues) const;

    void GetSolutionCoefficientVectorB(Vector& rValues) const;

    const GeometryType& GetGeometryMirror() const
    {
        return *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
    }
    
    array_1d<double, 3> mNormalParameterSpaceA = ZeroVector(3);
    array_1d<double, 3> mNormalParameterSpaceB = ZeroVector(3);

    array_1d<double, 3> mNormalPhysicalSpaceA = ZeroVector(3);
    array_1d<double, 3> mNormalPhysicalSpaceB = ZeroVector(3);
    ConstitutiveLaw::Pointer mpConstitutiveLaw;

    void ComputeTaylorExpansionContribution(const GeometryType& rGeometry, const Vector& rDistanceVector, Vector& H_sum_vec);
    void ComputeGradientTaylorExpansionContribution(const GeometryType& rGeometry, const Vector& rDistanceVector, Matrix& grad_H_sum);
    double ComputeTaylorTerm(double derivative, double dx, IndexType n_k, double dy, IndexType k);
    double ComputeTaylorTerm3D(double derivative, double dx, IndexType k_x, double dy, IndexType k_y, double dz, IndexType k_z);
    struct ConstitutiveVariables
    {
        ConstitutiveLaw::StrainVectorType StrainVector;
        ConstitutiveLaw::StressVectorType StressVector;
        ConstitutiveLaw::VoigtSizeMatrixType D;

        ConstitutiveVariables(const std::size_t StrainSize)
        {
            if (StrainVector.size() != StrainSize)
                StrainVector.resize(StrainSize);
            if (StressVector.size() != StrainSize)
                StressVector.resize(StrainSize);
            if (D.size1() != StrainSize || D.size2() != StrainSize)
                D.resize(StrainSize, StrainSize);
            noalias(StrainVector) = ZeroVector(StrainSize);
            noalias(StressVector) = ZeroVector(StrainSize);
            noalias(D)            = ZeroMatrix(StrainSize, StrainSize);
        }
    };

    void CalculateB(const GeometryType& rGeometry, Matrix& rB, Matrix& r_DN_DX) const;
    void ApplyConstitutiveLaw(
        std::size_t matSize,
        Vector& rStrain,
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rConstitutiVariables);

    unsigned int mDim = 0;
    IndexType mBasisFunctionsOrder = 0;
    bool mIsGapSbmCoupling = false;
    Vector mDistanceVectorB;
    double mPenalty = 0.0;
    double mNitschePenalty = 1.0;
};

} // namespace Kratos
