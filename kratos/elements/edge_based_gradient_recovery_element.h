//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#if !defined(KRATOS_EDGE_BASED_GRADIENT_RECOVIERY_ELEMENT_H_INCLUDED )
#define  KRATOS_EDGE_BASED_GRADIENT_RECOVIERY_ELEMENT_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"

namespace Kratos
{

///@addtogroup KratosCore
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// An element to compute the a nodal variable from an embedded skin
/**
 */

/**
 * @brief Element for the edge-based gradient recovery process
 * This element implements the edge-based gradient recovery technique
 * described in https://onlinelibrary.wiley.com/doi/epdf/10.1002/nme.4374.
 * It is intended to be used in combination with the KratosCore process
 * implemented in edge_based_gradient_recovery_process.h.
 * @tparam TDim Problem dimension
 */
template<std::size_t TDim>
class KRATOS_API(KRATOS_CORE) EdgeBasedGradientRecoveryElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EdgeBasedGradientRecoveryElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(EdgeBasedGradientRecoveryElement);

    /// Number of nodes of the element.
    constexpr static std::size_t NumNodes = 2;

    /// Single node block size
    constexpr static std::size_t BlockSize = TDim;

    /// Problem size
    constexpr static std::size_t LocalSize = NumNodes * BlockSize;

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

    /// Element types definition
    using IndexType = Element::IndexType;
    using DofsArrayType = Element::DofsArrayType;
    using DofsVectorType = Element::DofsVectorType;
    using EquationIdVectorType = Element::EquationIdVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    EdgeBasedGradientRecoveryElement(IndexType NewId = 0)
        : Element(NewId)
    {}

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    EdgeBasedGradientRecoveryElement(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes)
    {}

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    EdgeBasedGradientRecoveryElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {}

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    EdgeBasedGradientRecoveryElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~EdgeBasedGradientRecoveryElement() override = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
     * Returns a pointer to a new EdgeBasedGradientRecoveryElement element, created using given input
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<EdgeBasedGradientRecoveryElement<TDim>>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    /**
     * @brief Calculate the element's local contributions to the system for the current step.
     *
     * @param rLeftHandSideMatrix Reference to the local Left Hand Side (LHS) matrix
     * @param rRightHandSideVector Reference to the local Right Hand Side (RHS) vector
     * @param rCurrentProcessInfo Reference to the model part of interest ProcessInfo container
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Calculate the element's local Left Hand Side (LHS) contribution to the system for the current step.
     *
     * @param rLeftHandSideMatrix Reference to the local Left Hand Side (LHS) matrix
     * @param rCurrentProcessInfo Reference to the model part of interest ProcessInfo container
     */
    void CalculateLeftHandSide(
        MatrixType &rLeftHandSideMatrix,
        const ProcessInfo &rCurrentProcessInfo) override;

    /**
     * @brief Calculate the element's local Reft Hand Side (RHS) contribution to the system for the current step.
     *
     * @param rRightHandSideVector Reference to the local Right Hand Side (RHS) vector
     * @param rCurrentProcessInfo Reference to the model part of interest ProcessInfo container
     */
    void CalculateRightHandSide(
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override;

    /**
     * @brief Checks the input and that all required Kratos variables have been registered.
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @return 0 if no errors were found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     */
    int Check(const ProcessInfo &rCurrentProcessInfo) const override;

    /**
     * @brief Provides the global indices for each one of this element's local rows
     * This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * @brief Get the Dof List object
     * Returns a list of the element's Dofs
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Elemental Data
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "EdgeBasedGradientRecoveryElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "EdgeBasedGradientRecoveryElement";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {}

    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    EdgeBasedGradientRecoveryElement & operator=(EdgeBasedGradientRecoveryElement const& rOther) = delete;

    /// Copy constructor.
    EdgeBasedGradientRecoveryElement(EdgeBasedGradientRecoveryElement const& rOther) = delete;

    ///@}

}; // Class EdgeBasedGradientRecoveryElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<std::size_t TDim>
inline std::istream &operator>>(
    std::istream &rIStream,
    EdgeBasedGradientRecoveryElement<TDim> &rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TDim>
inline std::ostream &operator<<(
    std::ostream &rOStream,
    const EdgeBasedGradientRecoveryElement<TDim> &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // KratosCore group

} // namespace Kratos.

#endif // KRATOS_EDGE_BASED_GRADIENT_RECOVIERY_ELEMENT_H_INCLUDED  defined
