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

#if !defined(KRATOS_CONVECTION_DIFFUSION_REACTION_ELEMENT_H_INCLUDED)
#define KRATOS_CONVECTION_DIFFUSION_REACTION_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

template <
    unsigned int TDim,
    unsigned int TNumNodes,
    class TConvectionDiffusionReactionData>
class ConvectionDiffusionReactionElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = Element;

    /// Node type (default is: Node)
    using NodeType = Node;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

    /// Definition of nodes container type, redefined from GeometryType
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;

    /// Vector type for local contributions to the linear system
    using VectorType = Vector;

    /// Matrix type for local contributions to the linear system
    using MatrixType = Matrix;

    using IndexType = std::size_t;

    using EquationIdVectorType = std::vector<IndexType>;

    using DofsVectorType = std::vector<Dof<double>::Pointer>;

    /// Type for an array of shape function gradient matrices
    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    using PropertiesType = typename BaseType::PropertiesType;

    using ConvectionDiffusionReactionDataType = TConvectionDiffusionReactionData;


    using CurrentElementType =
        ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ConvectionDiffusionReactionElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ConvectionDiffusionReactionElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit ConvectionDiffusionReactionElement(
        IndexType NewId = 0)
    : Element(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    ConvectionDiffusionReactionElement(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    ConvectionDiffusionReactionElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    ConvectionDiffusionReactionElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        typename PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    ConvectionDiffusionReactionElement(
        ConvectionDiffusionReactionElement const& rOther)
    : Element(rOther)
    {
    }

    /**
     * Destructor
     */
    ~ConvectionDiffusionReactionElement() override = default;

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
        typename PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<CurrentElementType>(
            NewId, this->GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }

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
        typename PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<CurrentElementType>(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<CurrentElementType>(
            NewId, this->GetGeometry().Create(ThisNodes), this->pGetProperties());
        KRATOS_CATCH("");
    }

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& CurrentProcessInfo) const override;

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& CurrentProcessInfo) const override;

    void GetValuesVector(
        Vector& rValues,
        int Step = 0) const override;

    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0) const override;

    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0) const override;

    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide
     * methods they can be managed internally with a private method to do the
     * same calculations only once: MANDATORY
     */

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief CalculateLocalVelocityContribution Calculate the local contribution in terms of velocity and pressure.
     * @param rDampMatrix Local finite element system matrix (output)
     * @param rRightHandSideVector Local finite element residual vector (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    void CalculateLocalVelocityContribution(
        MatrixType& rDampingMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * ELEMENTS inherited from this class must implement this methods
     * if they need to add dynamic element contributions
     * note: second derivatives means the accelerations if the displacements are the dof of the analysis
     * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
     * CalculateSecondDerivativesContributions,
     * CalculateSecondDerivativesLHS, CalculateSecondDerivativesRHS methods are : OPTIONAL
     */

    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix: the elemental mass matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;
    /**
     * this is called during the assembling process in order
     * to calculate the elemental damping matrix
     * @param rDampingMatrix: the elemental damping matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     * this method is: MANDATORY
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ConvectionDiffusionReactionElement #" << Id();
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CDR" << TConvectionDiffusionReactionData::GetName();
    }

    ///@}

protected:
    ///@name Protected operations
    ///@{

    /**
     * @brief Get the Values Array
     *
     * @param rValues       Return values array
     * @param Step          Step
     */
    void GetValuesArray(
        BoundedVector<double, TNumNodes>& rValues,
        const int Step = 0) const;

    /**
     * @brief Calculates shape function data for this element
     *
     * @param rGaussWeights Gauss point weights list
     * @param rNContainer   Shape function values. Each row contains shape functions for respective gauss point
     * @param rDN_DX        List of matrices containing shape function derivatives for each gauss point
     */
    virtual void CalculateGeometryData(
        Vector& rGaussWeights,
        Matrix& rNContainer,
        ShapeFunctionDerivativesArrayType& rDN_DX) const;

    void AddLumpedMassMatrix(
        Matrix& rMassMatrix,
        const double Mass) const;

    void AddDampingMatrixGaussPointContributions(
        Matrix& rDampingMatrix,
        const double ReactionTerm,
        const double EffectiveKinematicViscosity,
        const Vector& rVelocityConvectiveTerms,
        const double GaussWeight,
        const Vector& rGaussShapeFunctions,
        const Matrix& rGaussdNa_dNb) const;

    ///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }

    ///@}

}; // Class ConvectionDiffusionReactionElement

///@}
///@name Input and output
///@{

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
inline std::istream& operator>>(
    std::istream& rIStream,
    ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>& rThis);

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_CONVECTION_DIFFUSION_REACTION_ELEMENT_H_INCLUDED defined
