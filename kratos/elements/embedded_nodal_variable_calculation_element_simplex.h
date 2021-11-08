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

#if !defined(KRATOS_EMBEDDED_NODAL_VARIABLE_CALCULATION_ELEMENT_H_INCLUDED )
#define  KRATOS_EMBEDDED_NODAL_VARIABLE_CALCULATION_ELEMENT_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"

// Application includes


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
template< class TVarType >
class KRATOS_API(KRATOS_CORE) EmbeddedNodalVariableCalculationElementSimplex : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EmbeddedNodalVariableCalculationElementSimplex
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(EmbeddedNodalVariableCalculationElementSimplex);

    /// Node type (default is: Node<3>)
    typedef Node <3> NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    /// Vector type for local contributions to the linear system
    typedef Vector VectorType;

    /// Matrix type for local contributions to the linear system
    typedef Matrix MatrixType;

    /// Element types definition
    typedef Element::IndexType IndexType;
    typedef Element::DofsArrayType DofsArrayType;
    typedef Element::DofsVectorType DofsVectorType;
    typedef Element::EquationIdVectorType EquationIdVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    EmbeddedNodalVariableCalculationElementSimplex(IndexType NewId = 0)
        : Element(NewId)
    {}

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    EmbeddedNodalVariableCalculationElementSimplex(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes)
    {}

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    EmbeddedNodalVariableCalculationElementSimplex(
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
    EmbeddedNodalVariableCalculationElementSimplex(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~EmbeddedNodalVariableCalculationElementSimplex() override = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
     * Returns a pointer to a new EmbeddedNodalVariableCalculationElementSimplex element, created using given input
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
        return Kratos::make_intrusive<EmbeddedNodalVariableCalculationElementSimplex<TVarType>>(NewId, GetGeometry().Create(ThisNodes), pProperties);
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
        buffer << "EmbeddedNodalVariableCalculationElementSimplex #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "EmbeddedNodalVariableCalculationElementSimplex";
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

    const array_1d<double, 2> GetDistanceBasedShapeFunctionValues();

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
    EmbeddedNodalVariableCalculationElementSimplex & operator=(EmbeddedNodalVariableCalculationElementSimplex const& rOther) = delete;

    /// Copy constructor.
    EmbeddedNodalVariableCalculationElementSimplex(EmbeddedNodalVariableCalculationElementSimplex const& rOther) = delete;

    ///@}

}; // Class EmbeddedNodalVariableCalculationElementSimplex

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template <class TVarType>
inline std::istream &operator>>(
    std::istream &rIStream,
    EmbeddedNodalVariableCalculationElementSimplex<TVarType> &rThis)
{
    return rIStream;
}

/// output stream function
template <class TVarType>
inline std::ostream &operator<<(
    std::ostream &rOStream,
    const EmbeddedNodalVariableCalculationElementSimplex<TVarType> &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // KratosCore group

} // namespace Kratos.

#endif // KRATOS_EMBEDDED_NODAL_VARIABLE_CALCULATION_ELEMENT_H_INCLUDED  defined
