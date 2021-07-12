//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application


#if !defined(KRATOS_TRUSS_ELEMENT_H_INCLUDED)
#define KRATOS_TRUSS_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"

// Application includes
#include "iga_application_variables.h"

namespace Kratos
{

class TrussElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( TrussElement );

    /// Size types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor using an array of nodes
    TrussElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {};

    /// Constructor using an array of nodes with properties
    TrussElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {};

    /// Default constructor necessary for serialization
    TrussElement()
        : Element()
    {};

    /// Destructor.
    ~TrussElement() override = default;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override
    {
        return Kratos::make_intrusive<TrussElement>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override
    {
        return Kratos::make_intrusive<TrussElement>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    ///@}
    ///@name Degrees of freedom
    ///@{

    /// Relates the degrees of freedom of the element geometry
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /// Sets the ID's of the element degrees of freedom
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    ///@}
    ///@name Analysis stages
    ///@{

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /// Computes RHS
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType nb_nodes = GetGeometry().size();
        const SizeType nb_dofs = nb_nodes * 3;

        if (rRightHandSideVector.size() != nb_dofs)
            rRightHandSideVector.resize(nb_dofs);
        noalias(rRightHandSideVector) = ZeroVector(nb_dofs);

        MatrixType left_hand_side_matrix;

        CalculateAll(left_hand_side_matrix, rRightHandSideVector,
            rCurrentProcessInfo, false, true);
    }

    /// Computes LHS
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType nb_nodes = GetGeometry().size();
        const SizeType nb_dofs = nb_nodes * 3;

        VectorType right_hand_side_vector;

        if (rLeftHandSideMatrix.size1() != nb_dofs && rLeftHandSideMatrix.size2() != nb_dofs)
            rLeftHandSideMatrix.resize(nb_dofs, nb_dofs);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(nb_dofs, nb_dofs);

        CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
            rCurrentProcessInfo, true, false);
    }

    /// Computes LHS and RHS
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType nb_nodes = GetGeometry().size();
        const SizeType nb_dofs = nb_nodes * 3;

        if (rRightHandSideVector.size() != nb_dofs)
            rRightHandSideVector.resize(nb_dofs);
        noalias(rRightHandSideVector) = ZeroVector(nb_dofs);

        if (rLeftHandSideMatrix.size1() != nb_dofs && rLeftHandSideMatrix.size2() != nb_dofs)
            rLeftHandSideMatrix.resize(nb_dofs, nb_dofs);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(nb_dofs, nb_dofs);

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
            rCurrentProcessInfo, true, true);
    }

    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool ComputeLeftHandSide,
        const bool ComputeRightHandSide);

    /// Updates the constitutive law
    void FinalizeSolutionStep(
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Load functions
    ///@{

    void CalculateBodyForces(Vector& rBodyForces);

    bool HasSelfWeight() const;

    ///@}
    ///@name Explicit dynamic functions
    ///@{

    void AddExplicitContribution(
        const VectorType& rRHSVector,
        const Variable<VectorType>& rRHSVariable,
        const Variable<double >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo) override;

    void AddExplicitContribution(
        const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
        const Variable<array_1d<double, 3>>& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Dynamic functions
    ///@{

    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Calculates the mass matrix with use of the lumped mass vector
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Calculates lumped mass vector
    void CalculateLumpedMassVector(
        VectorType &rLumpedMassVector,
        const ProcessInfo &rCurrentProcessInfo) const override;

    /// Get Displacemnts
    void GetValuesVector(
        Vector& rValues,
        int Step = 0) const override;

    /// Get Velocities
    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0) const override;

    /// Get Accelerations
    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0) const override;

    ///@}
    ///@name Output functions
    ///@{

    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Info
    ///@{

    /// Check provided parameters
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "\"TrussElement\" #" << Id()
            << " with geometry #" << this->GetGeometry().Id() << " with center in: "
            << this->GetGeometry().Center() << std::endl;
    }

    ///@}
private:
    ///@name Private member
    ///@{

    std::vector<array_1d<double, 3>> mReferenceBaseVector;

    /// The vector containing the constitutive laws for all integration points.
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    ///@}
    ///@name Private operations
    ///@{

    /// Initializes constitutive law vector and materials.
    void InitializeMaterial();

    /// Computes the base vector at a integration point position.
    array_1d<double, 3> CalculateActualBaseVector(IndexType IntegrationPointIndex) const;

    /// Computes Green Lagrange Strain for all integration points
    void CalculateGreenLagrangeStrain(
        std::vector<double>& rGreenLagrangeVector) const;

    void CalculateTangentModulus(
        std::vector<double>& rTangentModulusVector,
        const ProcessInfo& rCurrentProcessInfo);

    /// Computes prestress
    double CalculatePrestressPK2(double reference_a, double actual_a) const;

    /// Computes PK2 stress
    void CalculateStressPK2(
        std::vector<double>& rStressVector,
        const ProcessInfo& rCurrentProcessInfo) const;

    /// Computes cauchy stress
    void CalculateStressCauchy(
        std::vector<double>& rStressVector,
        const ProcessInfo& rCurrentProcessInfo) const;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TrussElement);
        rSerializer.save("ReferenceBaseVector", mReferenceBaseVector);
        rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TrussElement);
        rSerializer.load("ReferenceBaseVector", mReferenceBaseVector);
        rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
    }

    ///@}
};

} // namespace Kratos

#endif // !defined(KRATOS_TRUSS_ELEMENT_H_INCLUDED)
