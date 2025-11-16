//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_FLUID_ADJOINT_ELEMENT_H)
#define KRATOS_FLUID_ADJOINT_ELEMENT_H

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/constitutive_law.h"
#include "includes/element.h"
#include "includes/properties.h"
#include "utilities/adjoint_extensions.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
class FluidAdjointElement : public Element
{
    class ThisExtensions : public AdjointExtensions
    {
        Element* mpElement;

    public:
        explicit ThisExtensions(Element* pElement);

        void GetFirstDerivativesVector(
            std::size_t NodeId,
            std::vector<IndirectScalar<double>>& rVector,
            std::size_t Step) override;

        void GetSecondDerivativesVector(
            std::size_t NodeId,
            std::vector<IndirectScalar<double>>& rVector,
            std::size_t Step) override;

        void GetAuxiliaryVector(
            std::size_t NodeId,
            std::vector<IndirectScalar<double>>& rVector,
            std::size_t Step) override;

        void GetFirstDerivativesVariables(std::vector<VariableData const*>& rVariables) const override;

        void GetSecondDerivativesVariables(std::vector<VariableData const*>& rVariables) const override;

        void GetAuxiliaryVariables(std::vector<VariableData const*>& rVariables) const override;
    };

public:
    ///@name Type Definitions
    ///@{

    using BaseType = Element;

    using NodesArrayType = typename BaseType::NodesArrayType;

    using PropertiesType = typename BaseType::PropertiesType;

    using GeometryType = typename BaseType::GeometryType;

    using VectorType = typename BaseType::VectorType;

    using MatrixType = typename BaseType::MatrixType;

    using IndexType = std::size_t;

    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    constexpr static IndexType TBlockSize = TDim + 1;

    constexpr static IndexType TElementLocalSize = TBlockSize * TNumNodes;

    constexpr static IndexType TCoordsLocalSize = TDim * TNumNodes;

    using VectorF = BoundedVector<double, TElementLocalSize>;

    KRATOS_CLASS_POINTER_DEFINITION(FluidAdjointElement);

    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    FluidAdjointElement(IndexType NewId = 0);

    /**
     * Constructor using Geometry
     */
    FluidAdjointElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    /**
     * Constructor using Properties
     */
    FluidAdjointElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);

    /**
     * Destructor
     */
    ~FluidAdjointElement() override;

    ///@}
    ///@name Operations
    ///@{

    /**
     * ELEMENTS inherited from this class have to implement next
     * Create and Clone methods: MANDATORY
     */

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override;

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override;

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& ThisNodes) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult the elemental equation ID vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rElementalEquationIdList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /// Returns the adjoint values stored in this element's nodes.
    void GetValuesVector(
        VectorType& rValues,
        int Step = 0) const override;

    /// Returns the adjoint velocity values stored in this element's nodes.
    void GetFirstDerivativesVector(
        VectorType& rValues,
        int Step = 0) const override;

    void GetSecondDerivativesVector(
        VectorType& rValues,
        int Step) const override;

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
        MatrixType &rDampMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override;

    void CalculateFirstDerivativesLHS(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateSecondDerivativesLHS(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateSensitivityMatrix(
        const Variable<array_1d<double, 3>>& rSensitivityVariable,
        Matrix& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(
        const Variable<Vector>& rVariable,
        Vector& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(
        const Variable<Matrix>& rVariable,
        Matrix& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

protected:
    ///@name Protected Members
    ///@{

    ConstitutiveLaw::Pointer mpConstitutiveLaw = nullptr;

    ///@}
    ///@name Protected Operations
    ///@{

    void AddFluidResidualsContributions(
        VectorType& rResidual,
        const ProcessInfo& rCurrentProcessInfo);

    void AddFluidFirstDerivatives(
        MatrixType& rDerivativesMatrix,
        const ProcessInfo& rCurrentProcessInfo,
        const double MassTermsDerivativesWeight = 1.0);

    void AddFluidSecondDerivatives(
        MatrixType& rDerivativesMatrix,
        const ProcessInfo& rCurrentProcessInfo);

    void AddFluidShapeDerivatives(
        MatrixType& rDerivativesMatrix,
        const ProcessInfo& rCurrentProcessInfo);

    void CalculateGeometryData(
        Vector& rGaussWeights,
        Matrix& rNContainer,
        ShapeFunctionDerivativesArrayType& rDN_DX,
        const GeometryData::IntegrationMethod& rIntegrationMethod) const;

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_ADJOINT_ELEMENT_H