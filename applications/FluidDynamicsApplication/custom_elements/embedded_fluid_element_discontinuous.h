//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//  Co-authors:      Jordi Cotela
//

#ifndef KRATOS_EMBEDDED_FLUID_ELEMENT_DISCONTINUOUS_H
#define KRATOS_EMBEDDED_FLUID_ELEMENT_DISCONTINUOUS_H

#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "modified_shape_functions/modified_shape_functions.h"

#include "includes/cfd_variables.h"
#include "custom_elements/fluid_element.h"

#include "custom_utilities/embedded_discontinuous_data.h"

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
class EmbeddedFluidElementDiscontinuous : public TBaseElement
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EmbeddedFluidElementDiscontinuous
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(EmbeddedFluidElementDiscontinuous);

    /// Node type (default is: Node<3>)
    typedef Node<3> NodeType;

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

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    /// Type for shape function values container
    typedef Kratos::Vector ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef Kratos::Matrix ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef Geometry<NodeType>::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    constexpr static unsigned int Dim = TBaseElement::Dim;
    constexpr static unsigned int NumNodes = TBaseElement::NumNodes;
    constexpr static unsigned int BlockSize = TBaseElement::BlockSize;
    constexpr static unsigned int LocalSize = TBaseElement::LocalSize;
    constexpr static unsigned int StrainSize = TBaseElement::StrainSize;

    using BaseElementData = typename TBaseElement::ElementData;
    using EmbeddedDiscontinuousElementData = EmbeddedDiscontinuousData< BaseElementData >;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    EmbeddedFluidElementDiscontinuous(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    EmbeddedFluidElementDiscontinuous(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    EmbeddedFluidElementDiscontinuous(IndexType NewId, Geometry<NodeType>::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    EmbeddedFluidElementDiscontinuous(IndexType NewId, Geometry<NodeType>::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    ~EmbeddedFluidElementDiscontinuous() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Create a new element of this type
    /**
     * Returns a pointer to a new EmbeddedFluidElementDiscontinuous element, created using given input
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
    /** For EmbeddedFluidElementDiscontinuous, this initializes the discontinuous
     * level set (ELEMENTAL_DISTANCES) and the nodal imposed velocity (EMBEDDED_VELOCITY)
     */
    void Initialize(const ProcessInfo &rCurrentProcessInfo) override;

    /// Calculates both LHS and RHS contributions
    /**
     * Computes the LHS and RHS elementar matrices. If the element is split
     * includes the contribution of the level set boundary condition imposition.
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rRightHandSideVector reference to the RHS vector
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
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
    void InitializeGeometryData(EmbeddedDiscontinuousElementData& rData) const;

    /**
     * @brief Non-intersected element geometry data fill
     * This method sets the data structure geometry fields (shape functions, gradients, ...) for a non-intersected element.
     * @param rData reference to the element data structure
     */
    void DefineStandardGeometryData(EmbeddedDiscontinuousElementData& rData) const;

    /**
     * @brief Intersected element geometry data fill
     * This method sets the data structure geometry fields (shape functions, gradients, interface normals, ...) for an
     * intersected element. To do that, the modified shape functions utility is firstly created and then called to perform
     * all operations in both the positive and negative sides of the element.
     * @param rData reference to the element data structure
     */
    void DefineCutGeometryData(EmbeddedDiscontinuousElementData& rData) const;

    /**
     * @brief For an intersected element, normalize the interface normals
     * This method normalizes the interface normals for an intersected element.
     * @param rNormals interface normals container
     * @param Tolerance tolerance to avoid division by 0 when normalizing
     */
    void NormalizeInterfaceNormals(
        typename EmbeddedDiscontinuousElementData::InterfaceNormalsType& rNormals,
        double Tolerance) const;

    /**
    * This method adds the no-penetration condition penalty level set contribution.
    * @param rLHS reference to the LHS matrix
    * @param rRHS reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddNormalPenaltyContribution(
        MatrixType& rLHS,
        VectorType& rRHS,
        const EmbeddedDiscontinuousElementData& rData) const;

    /**
    * This method adds the no-penetration condition adjoint term level set contribution.
    * @param rLHS reference to the LHS matrix
    * @param rRHS reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddNormalSymmetricCounterpartContribution(
        MatrixType& rLHS,
        VectorType& rRHS,
        const EmbeddedDiscontinuousElementData& rData) const;

    /**
    * This method adds the tangential stress condition penalty level set contribution.
    * @param rLHS reference to the LHS matrix
    * @param rRHS reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddTangentialPenaltyContribution(
        MatrixType& rLHS,
        VectorType& rRHS,
        const EmbeddedDiscontinuousElementData& rData) const;

    /**
    * This method adds the tangential stress condition adjoint term level set contribution.
    * @param rLHS reference to the LHS matrix
    * @param rRHS reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddTangentialSymmetricCounterpartContribution(
        MatrixType& rLHS,
        VectorType& rRHS,
        const EmbeddedDiscontinuousElementData& rData) const;

    /**
     * This method computes the penalty coefficient for the Nitsche normal imposition
     * @param rData reference to element data structure
     * @param rN the current Gauss pt. shape functions vector
     * @return double The normal penalty coefficient value
     */
    double ComputeNormalPenaltyCoefficient(
        const EmbeddedDiscontinuousElementData& rData,
        const Vector& rN) const;

    /**
     * This method computes the Nitsche coefficients for the Nitsche normal imposition
     * @param rData reference to element data structure
     * @return a pair of double containing the two coefficients
     */
    std::pair<const double, const double> ComputeTangentialPenaltyCoefficients(const EmbeddedDiscontinuousElementData& rData) const;

    /**
     * This method computes the Nitsche coefficients for the Nitsche tangential imposition
     * @param rData reference to element data structure
     * @return a pair of double containing the two coefficients
     */
    std::pair<const double, const double> ComputeTangentialNitscheCoefficients(const EmbeddedDiscontinuousElementData& rData) const;


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
     * @brief Calculates the drag force
     * For an intersected element, this method calculates the drag force.
     * Note that the drag force includes both the shear and the pressure contributions.
     * @param rData reference to the embedded elemental data
     * @param rDragForce reference to the computed drag force
     */
    void CalculateDragForce(
        EmbeddedDiscontinuousElementData& rData,
        array_1d<double,3>& rDragForce) const;

    /**
     * @brief Calculates the location of the drag force
     * For an intersected element, this method calculates the drag force location.
     * Note that the drag force includes both the shear and the pressure contributions.
     * @param rData reference to the embedded elemental data
     * @param rDragForce reference to the computed drag force
     */
    void CalculateDragForceCenter(
        EmbeddedDiscontinuousElementData& rData,
        array_1d<double,3>& rDragForceLocation) const;

    /**
     * @brief Auxiliary method to get the density value
     * This auxiliary method interfaces the density get in order to make possible the
     * use of the embedded element with both property-based and nodal-based density formulations.
     * For the standard case (property-based formulations) the method is not specialized.
     * In case a nodal density base formulation is used, it needs to be specialized.
     * @param rData Embedded element data container
     * @param NodeIndex The local index node for which the density is retrieved (only used in nodal density formulations)
     * @return double The density value
     */
    inline double AuxiliaryDensityGetter(
        const EmbeddedDiscontinuousElementData& rData,
        const unsigned int NodeIndex) const;

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
    EmbeddedFluidElementDiscontinuous& operator=(EmbeddedFluidElementDiscontinuous const& rOther);

    /// Copy constructor.
    EmbeddedFluidElementDiscontinuous(EmbeddedFluidElementDiscontinuous const& rOther);

    ///@}


}; // Class EmbeddedFluidElementDiscontinuous

namespace EmbeddedDiscontinuousInternals {

template <size_t TDim, size_t TNumNodes>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator(
    const Element &rElement,
    const Vector &rElementalDistances);

template <size_t TDim, size_t TNumNodes>
ModifiedShapeFunctions::Pointer GetContinuousShapeFunctionCalculator(
    const Element &rElement,
    const Vector &rElementalDistances);
}

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
    EmbeddedFluidElementDiscontinuous<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(
    std::ostream& rOStream,
    const EmbeddedFluidElementDiscontinuous<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_EMBEDDED_FLUID_ELEMENT_DISCONTINUOUS_H
