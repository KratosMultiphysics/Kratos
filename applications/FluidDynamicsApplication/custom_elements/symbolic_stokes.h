//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:  Mohit Tyagi, Daniel Diez, Pablo Becker
//  Co-authors:
//

#if !defined(KRATOS_SYMBOLIC_STOKES)
#define  KRATOS_SYMBOLIC_STOKES

// System includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "includes/cfd_variables.h"
#include "custom_elements/fluid_element.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

/*The "SymbolicStokes" element is an element based on the Variational Multiscale Stabilization technique (VMS).
*/

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

template< class TElementData >
class SymbolicStokes : public FluidElement<TElementData>
{
public:

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SymbolicStokes);

    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef std::vector<std::size_t> EquationIdVectorType;
    typedef std::vector< Dof<double>::Pointer > DofsVectorType;
    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;
    typedef typename FluidElement<TElementData>::ShapeFunctionsType ShapeFunctionsType;
    typedef typename FluidElement<TElementData>::ShapeFunctionDerivativesType ShapeFunctionDerivativesType;
    typedef typename FluidElement<TElementData>::ShapeFunctionDerivativesArrayType ShapeFunctionDerivativesArrayType;
    constexpr static unsigned int Dim = FluidElement<TElementData>::Dim;
    constexpr static unsigned int NumNodes = FluidElement<TElementData>::NumNodes;
    constexpr static unsigned int BlockSize = FluidElement<TElementData>::BlockSize;
    constexpr static unsigned int LocalSize = FluidElement<TElementData>::LocalSize;
    constexpr static unsigned int StrainSize = (Dim - 1) * 3;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constuctor.
    /**
    * @param NewId Index number of the new element (optional)
    */
    SymbolicStokes(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
    * @param NewId Index of the new element
    * @param ThisNodes An array containing the nodes of the new element
    */
    SymbolicStokes(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
    * @param NewId Index of the new element
    * @param pGeometry Pointer to a geometry object
    */
    SymbolicStokes(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
    * @param NewId Index of the new element
    * @param pGeometry Pointer to a geometry object
    * @param pProperties Pointer to the element's properties
    */
    SymbolicStokes(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    virtual ~SymbolicStokes();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
    * Returns a pointer to a new SymbolicStokes element, created using given input.
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
    * Returns a pointer to a new FluidElement element, created using given input.
    * @param NewId the ID of the new element
    * @param pGeom a pointer to the geomerty to be used to create the element
    * @param pProperties the properties assigned to the new element
    * @return a Pointer to the new element
    */
    Element::Pointer Create(IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Function to visualize the divergence field

    /**
     * @brief Get the Value On Integration Points object (used to visualize the divergence field)
     *
     * @param rVariable Variable to be retrieved (implementation supports DIVERGENCE)
     * @param rValues Vector for the values at the Gauss integration points
     * @param rCurrentProcessInfo ProcessInfo object
     */
    void GetValueOnIntegrationPoints(
        const Variable<double> &rVariable,
        std::vector<double> &rValues,
        const ProcessInfo &rCurrentProcessInfo) override;


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
     * @brief Computes time integrated LHS and RHS arrays
     * This method computes both the Left Hand Side and
     * Right Hand Side time integrated contributions.
     * @param rData Reference to the element data container
     * @param rLHS Reference to the Left Hand Side matrix to be filled
     * @param rRHS Reference to the Right Hand Side vector to be filled
     */
    void AddTimeIntegratedSystem(
        TElementData& rData,
        MatrixType& rLHS,
        VectorType& rRHS) override;

    /**
     * @brief Computes the time integrated LHS matrix
     * This method computes the Left Hand Side time integrated contribution
     * @param rData Reference to the element data container
     * @param rLHS Reference to the Left Hand Side matrix to be filled
     */
    void AddTimeIntegratedLHS(
        TElementData& rData,
        MatrixType& rLHS) override;

    /**
     * @brief Computes the time integrated RHS vector
     * This method computes the Right Hand Side time integrated contribution
     * @param rData Reference to the element data container
     * @param rRHS Reference to the Right Hand Side matrix to be filled
     */
    void AddTimeIntegratedRHS(
        TElementData& rData,
        VectorType& rRHS) override;

    /**
     * @brief Computes the LHS Gauss pt. contribution
     * This method computes the contribution to the LHS of a Gauss pt.
     * @param rData Reference to the element data container
     * @param rLHS Reference to the Left Hand Side matrix to be filled
     */
    void ComputeGaussPointLHSContribution(
        TElementData& rData,
        MatrixType& rLHS);

    /**
     * @brief Computes the RHS Gaus  pt. contribution
     * This method computes the contribution to the RHS of a Gauss pt.
     * @param rData Reference to the element data container
     * @param rRHS Reference to the Right Hand Side vector to be filled
     */
    void ComputeGaussPointRHSContribution(
        TElementData& rData,
        VectorType& rRHS);

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
private:
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
    SymbolicStokes& operator=(SymbolicStokes const& rOther);

    /// Copy constructor.
    SymbolicStokes(SymbolicStokes const& rOther);

    ///@}

}; // Class SymbolicStokes
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template< class TElementData >
inline std::istream& operator >> (std::istream& rIStream,
    SymbolicStokes<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(std::ostream& rOStream,
    const SymbolicStokes<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_SYMBOLIC_STOKES
