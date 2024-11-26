//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Franziska Wahl
//

#pragma once

#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "modified_shape_functions/modified_shape_functions.h"

#include "includes/cfd_variables.h"
#include "custom_elements/fluid_element.h"

#include "custom_utilities/embedded_data.h"

namespace Kratos
{

///@addtogroup FluidDynamicsApplication
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

template< class TBaseElement >
class ShiftedBoundaryFluidElement : public TBaseElement
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ShiftedBoundaryFluidElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ShiftedBoundaryFluidElement);

    /// Node type (default is: Node)
    typedef Node NodeType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    /// Vector type for local contributions to the linear system
    typedef Vector VectorType;

    /// Matrix type for local contributions to the linear system
    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet< Dof<double>, IndexedObject > DofsArrayType;

    constexpr static std::size_t Dim = TBaseElement::Dim;
    constexpr static std::size_t NumNodes = TBaseElement::NumNodes;
    constexpr static std::size_t BlockSize = TBaseElement::BlockSize;
    constexpr static std::size_t LocalSize = TBaseElement::LocalSize;
    constexpr static std::size_t StrainSize = TBaseElement::StrainSize;

    using BaseElementData = typename TBaseElement::ElementData;
    using ShiftedBoundaryElementData = BaseElementData;  //NOTE EmbeddedData<BaseElementData> for (dis-)continuous level set

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    ShiftedBoundaryFluidElement(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    ShiftedBoundaryFluidElement(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    ShiftedBoundaryFluidElement(IndexType NewId, Geometry<NodeType>::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    ShiftedBoundaryFluidElement(IndexType NewId, Geometry<NodeType>::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    ~ShiftedBoundaryFluidElement() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Create a new element of this type
    /**
     * Returns a pointer to a new ShiftedBoundaryFluidElement element, created using given input
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            Properties::Pointer pProperties) const override;

    /// Create a new element of this type using given geometry
    /**
     * Returns a pointer to a new FluidElement element, created using given input
     * @param NewId the ID of the new element
     * @param pGeom a pointer to the geomerty to be used to create the element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            Geometry<NodeType>::Pointer pGeom,
                            Properties::Pointer pProperties) const override;

    /// Set up the element for solution.
    /** For discontinuous level set instead of point-based, this needs to initialize the discontinuous
     * level set distances (ELEMENTAL_DISTANCES) and the nodal imposed velocity (EMBEDDED_VELOCITY)
     */
    void Initialize(const ProcessInfo &rCurrentProcessInfo) override;

    /**
     * Computes the LHS and RHS elemental matrices.
     * If the element is flagged as INTERFACE it contains surrogate boundary faces, for which the boundary flux contribution is added.
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rRightHandSideVector reference to the RHS vector
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Computes the LHS elemental matrices.
     * If the element is flagged as INTERFACE it contains surrogate boundary faces, for which the boundary flux contribution is added.
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Computes the RHS elemental vector.
     * If the element is flagged as INTERFACE it contains surrogate boundary faces, for which the boundary flux contribution is added.
     * @param rRightHandSideVector reference to the RHS vector
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Computes an elemental double value
    /**
     * Given a double variable, this function computes its value inside de element.
     * If the function has not implemented this variable computation, calls the base class one.
     * @param rVariable Variable to be computed
     * @param rOutput Reference to the output double
     * @param rCurrentProcessInfo Reference to the process info
     */
    void Calculate(
        const Variable<double> &rVariable,
        double &rOutput,
        const ProcessInfo &rCurrentProcessInfo) override;

    /// Computes an elemental 3 components array value
    /**
     * Given a 3 components array variable, this function computes its value inside de element.
     * If the function has not implemented this variable computation, calls the base class one.
     * @param rVariable Variable to be computed
     * @param rOutput Reference to the output array
     * @param rCurrentProcessInfo Reference to the process info
     */
    void Calculate(
        const Variable<array_1d<double, 3>> &rVariable,
        array_1d<double, 3> &rOutput,
        const ProcessInfo &rCurrentProcessInfo) override;

    /// Computes an elemental vector value
    /**
     * Given a vector variable, this function computes its value inside de element.
     * If the function has not implemented this variable computation, calls the base class one.
     * @param rVariable Variable to be computed
     * @param rOutput Reference to the output vector
     * @param rCurrentProcessInfo Reference to the process info
     */
    void Calculate(
        const Variable<Vector> &rVariable,
        Vector &rOutput,
        const ProcessInfo &rCurrentProcessInfo) override;

    /// Computes an elemental matrix value
    /**
     * Given a matrix variable, this function computes its value inside de element.
     * If the function has not implemented this variable computation, calls the base class one.
     * @param rVariable Variable to be computed
     * @param rOutput Reference to the output matrix
     * @param rCurrentProcessInfo Reference to the process info
     */
    void Calculate(
        const Variable<Matrix> &rVariable,
        Matrix &rOutput,
        const ProcessInfo &rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    int Check(const ProcessInfo &rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    const Parameters GetSpecifications() const override;

    /// Turn back information as a string.
    std::string Info() const override;


    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;


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

    /**
     * @brief Current element data structure initialization
     * This method checks if the element is intersected and calls the elemental data filling methods accordingly.
     * @param rData reference to the element data structure
     */
    void InitializeGeometryData(ShiftedBoundaryElementData& rData) const;

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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Get the Surrogate Faces local Ids
     * This method returns a list with the local ids of the surrogate faces of INTERFACE elements.
     * Surrogate faces are faces shared with BOUNDARY elements and for which all their nodes  lie in the surrogate boundary.
     * @return std::vector<std::size_t> List with the surrogate faces local ids
     */
    std::vector<std::size_t> GetSurrogateFacesIds();

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
    ShiftedBoundaryFluidElement& operator=(ShiftedBoundaryFluidElement const& rOther);

    /// Copy constructor.
    ShiftedBoundaryFluidElement(ShiftedBoundaryFluidElement const& rOther);

    ///@}


}; // Class ShiftedBoundaryFluidElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TElementData >
inline std::istream& operator >>(
    std::istream& rIStream,
    ShiftedBoundaryFluidElement<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(
    std::ostream& rOStream,
    const ShiftedBoundaryFluidElement<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.
