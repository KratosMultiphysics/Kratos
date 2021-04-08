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

#if !defined(KRATOS_CONVECTION_DIFFUSION_REACTION_CROSS_WIND_STABILIZED_ELEMENT_H_INCLUDED)
#define KRATOS_CONVECTION_DIFFUSION_REACTION_CROSS_WIND_STABILIZED_ELEMENT_H_INCLUDED

// System includes
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"

// Application includes
#include "custom_elements/convection_diffusion_reaction_element.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim,
          unsigned int TNumNodes,
          class TConvectionDiffusionReactionData>
class ConvectionDiffusionReactionCrossWindStabilizedElement
: public ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>;

    /// Node type (default is: Node<3>)
    using NodeType = Node<3>;

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

    using PropertiesType = typename BaseType::PropertiesType;

    /// Type for an array of shape function gradient matrices
    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    using ConvectionDiffusionReactionDataType = TConvectionDiffusionReactionData;

    using CurrentElementType =
        ConvectionDiffusionReactionCrossWindStabilizedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ConvectionDiffusionReactionCrossWindStabilizedElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ConvectionDiffusionReactionCrossWindStabilizedElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit ConvectionDiffusionReactionCrossWindStabilizedElement(
        IndexType NewId = 0)
    : BaseType(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    ConvectionDiffusionReactionCrossWindStabilizedElement(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
    : BaseType(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    ConvectionDiffusionReactionCrossWindStabilizedElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    ConvectionDiffusionReactionCrossWindStabilizedElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    ConvectionDiffusionReactionCrossWindStabilizedElement(
        ConvectionDiffusionReactionCrossWindStabilizedElement const& rOther)
    : BaseType(rOther)
    {
    }

    /**
     * Destructor
     */
    ~ConvectionDiffusionReactionCrossWindStabilizedElement() override = default;

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
            NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
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
            NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
        KRATOS_CATCH("");
    }

    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide
     * methods they can be managed internally with a private method to do the
     * same calculations only once: MANDATORY
     */

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

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ConvectionDiffusionReactionCrossWindStabilizedElement #" << this->Id();
        return buffer.str();
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CDRCrossWind" << TConvectionDiffusionReactionData::GetName();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    /**
     * @brief Calculate gradient matrix for a vector
     *
     * Calculates the gradient matrix for a given vector variable.
     *
     * @param rOutput            Output matrix, rows contain the given vector indices, columns containt physical coordinate dimensions
     * @param rVariable          Vector variable
     * @param rShapeDerivatives  Shape function derivatives at the gauss point
     * @param Step               Time step
     */
    void CalculateGradient(
        BoundedMatrix<double, TDim, TDim>& rOutput,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rShapeDerivatives,
        const int Step = 0) const;

    /**
     * @brief Calculate gradient vector for a scalar
     *
     * Calculates the gradient vector for a given scalar variable.
     *
     * @param rOutput            Output vector
     * @param rVariable          Scalar variable
     * @param rShapeDerivatives  Shape function derivatives at the gauss point
     * @param Step               Time step
     */
    void CalculateGradient(
        array_1d<double, 3>& rOutput,
        const Variable<double>& rVariable,
        const Matrix& rShapeDerivatives,
        const int Step = 0) const;

    void CalculateContravariantMetricTensor(
        BoundedMatrix<double, TDim, TDim>& rOutput,
        const Matrix& rParameterDerivatives) const;

    virtual double GetDeltaTime(
        const ProcessInfo& rProcessInfo) const;

    double GetScalarVariableGradientNorm(
        const Matrix& rShapeFunctionDerivatives,
        const int Step = 0) const;

    /**
     * @brief Get the Geometry Parameter Derivatives object
     *
     * This method calculates partial derivatives of parametric coordinates(i.e. $\underline{\xi}$) of element
     * w.r.t. physical coordinates (i.e. $\underline{x}$)
     *
     * \[
     *      \frac{\partial \underline{\xi}}{\partial \underline{x}}
     * \]
     *
     * @return ShapeFunctionDerivativesArrayType
     */
    ShapeFunctionDerivativesArrayType GetGeometryParameterDerivatives() const;

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
}; // Class ConvectionDiffusionReactionCrossWindStabilizedElement

///@}
///@name Input and output
///@{

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
inline std::istream& operator>>(
    std::istream& rIStream,
    ConvectionDiffusionReactionCrossWindStabilizedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>& rThis);

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ConvectionDiffusionReactionCrossWindStabilizedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_CONVECTION_DIFFUSION_REACTION_CROSS_WIND_STABILIZED_ELEMENT_H_INCLUDED defined
