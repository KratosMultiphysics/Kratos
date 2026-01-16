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
#include "iga_application_variables.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

/**
 * @brief Interface coupling condition for fluid fields.
 *
 * Weakly enforces velocity continuity across adjacent patches using a
 * Nitsche/penalty-style coupling on the shared interface. Pressure DOFs
 * are accounted for in the DOF list and equation IDs but are not directly
 * coupled in this condition.
 */
class KRATOS_API(IGA_APPLICATION) FluidCouplingCondition : public Condition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(FluidCouplingCondition);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// Type for shape function derivatives container
    typedef Kratos::Matrix ShapeDerivativesType;

    FluidCouplingCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    FluidCouplingCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    FluidCouplingCondition() = default;
    ~FluidCouplingCondition() override = default;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<FluidCouplingCondition>(NewId, pGeom, pProperties);
    }

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<FluidCouplingCondition>(
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
        buffer << "\"FluidCouplingCondition\" #" << Id();
        return buffer.str();
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"FluidCouplingCondition\" #" << Id();
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

    void CalculateB(
        Matrix& rB,
        const ShapeDerivativesType& r_DN_DX) const;

    void GetVelocityCoefficientVectorA(Vector& rValues) const;
    void GetVelocityCoefficientVectorB(Vector& rValues) const;

private:

    const GeometryType& GetGeometryMirror() const
    {
        KRATOS_ERROR_IF_NOT(this->Has(NEIGHBOUR_GEOMETRIES))
            << "FluidCouplingCondition #" << this->Id() << " missing NEIGHBOUR_GEOMETRIES." << std::endl;
        const auto& neigh = this->GetValue(NEIGHBOUR_GEOMETRIES);
        KRATOS_ERROR_IF(neigh.empty())
            << "FluidCouplingCondition #" << this->Id() << ": NEIGHBOUR_GEOMETRIES is empty." << std::endl;
        return *neigh[0];
    }
    
    // Small helper to hold CL variables (2D Voigt size = 3)
    struct ConstitutiveVariables
    {
        ConstitutiveLaw::StrainVectorType StrainVector;
        ConstitutiveLaw::StressVectorType StressVector;
        ConstitutiveLaw::VoigtSizeMatrixType D;

        explicit ConstitutiveVariables(const SizeType StrainSize)
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

    // Applies the constitutive law at the interface using the given B-matrix
    // to compute strains. Fills stress vector and constitutive tensor in
    // rConstitutiveVariables via rValues and the stored constitutive law.
    // Apply CL using a specific geometry (A or B) so nodal values and sizes
    // are taken from the correct side of the interface
    void ApplyConstitutiveLaw(
        const GeometryType& rGeometry,
        const Matrix& rB,
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rConstitutiveVariables) const;

    // Constitutive law pointer (cloned from Properties)
    ConstitutiveLaw::Pointer mpConstitutiveLaw;

    array_1d<double, 3> mNormalParameterSpaceA = ZeroVector(3);
    array_1d<double, 3> mNormalParameterSpaceB = ZeroVector(3);

    array_1d<double, 3> mNormalPhysicalSpaceA = ZeroVector(3);
    array_1d<double, 3> mNormalPhysicalSpaceB = ZeroVector(3);

    // Polynomial order used to scale penalty (p^2/h)
    IndexType mBasisFunctionsOrder = 1;

    // Spatial velocity dimension inferred from geometry (2 or 3)
    SizeType mDim = 2;
    
    double mPenalty;
    double mNitschePenalty;
};

} // namespace Kratos
