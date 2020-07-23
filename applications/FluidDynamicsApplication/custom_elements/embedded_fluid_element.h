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

#ifndef KRATOS_EMBEDDED_FLUID_ELEMENT_H
#define KRATOS_EMBEDDED_FLUID_ELEMENT_H

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
class EmbeddedFluidElement : public TBaseElement
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EmbeddedFluidElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(EmbeddedFluidElement);

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
    using EmbeddedElementData = EmbeddedData< BaseElementData >;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    EmbeddedFluidElement(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    EmbeddedFluidElement(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    EmbeddedFluidElement(IndexType NewId, Geometry<NodeType>::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    EmbeddedFluidElement(IndexType NewId, Geometry<NodeType>::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    ~EmbeddedFluidElement() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Create a new element of this type
    /**
     * Returns a pointer to a new EmbeddedFluidElement element, created using given input
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
    /** For EmbeddedFluidElement, this initializes the nodal imposed velocity (EMBEDDED_VELOCITY)
     */
    void Initialize() override;

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
        ProcessInfo& rCurrentProcessInfo) override;

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

    /**
     * @brief Base element GetValueOnIntegrationPoints
     * Called to avoid reimplementing the variable types specializations not required
     */
    using TBaseElement::GetValueOnIntegrationPoints;

    /**
     * @brief Get the Value On Integration Points object
     * Computes the value in the Gauss pts. for a three component array variable
     * @param rVariable Array variable to be computed
     * @param rValues Computed gauss point values
     * @param rCurrentProcessInfo Current process info
     */
    void GetValueOnIntegrationPoints(
        const Variable<array_1d<double, 3>> &rVariable,
        std::vector<array_1d<double, 3>> &rValues,
        const ProcessInfo &rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{

    int Check(const ProcessInfo &rCurrentProcessInfo) override;

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

    void InitializeGeometryData(EmbeddedElementData& rData) const;

    void DefineStandardGeometryData(EmbeddedElementData& rData) const;

    void DefineCutGeometryData(EmbeddedElementData& rData) const;

    void NormalizeInterfaceNormals(typename EmbeddedElementData::InterfaceNormalsType& rNormals, double Tolerance) const;

    /**
    * This functions adds the no-penetration condition penalty level set contribution.
    * @param rLHS reference to the LHS matrix
    * @param rRHS reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddSlipNormalPenaltyContribution(
        MatrixType& rLHS,
        VectorType& rRHS,
        const EmbeddedElementData& rData) const;

    /**
    * This functions adds the no-penetration condition adjoint term level set contribution.
    * @param rLHS reference to the LHS matrix
    * @param rRHS reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddSlipNormalSymmetricCounterpartContribution(
        MatrixType& rLHS,
        VectorType& rRHS,
        const EmbeddedElementData& rData) const;

    /**
    * This functions adds the tangential stress condition penalty level set contribution.
    * @param rLHS reference to the LHS matrix
    * @param rRHS reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddSlipTangentialPenaltyContribution(
        MatrixType& rLHS,
        VectorType& rRHS,
        const EmbeddedElementData& rData) const;

    /**
    * This functions adds the tangential stress condition adjoint term level set contribution.
    * @param rLHS reference to the LHS matrix
    * @param rRHS reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddSlipTangentialSymmetricCounterpartContribution(
        MatrixType& rLHS,
        VectorType& rRHS,
        const EmbeddedElementData& rData) const;

    /**
     * This function computes the penalty coefficient for the Nitsche normal imposition
     * @param rData reference to element data structure
     */
    double ComputeSlipNormalPenaltyCoefficient(
        const EmbeddedElementData& rData) const;

    /**
     * This function computes the Nitsche coefficients for the Nitsche normal imposition
     * @param rData reference to element data structure
     * @return a pair of double containing the two coefficients
     */
    std::pair<const double, const double> ComputeSlipTangentialPenaltyCoefficients(
        const EmbeddedElementData& rData) const;

    /**
     * This function computes the Nitsche coefficients for the Nitsche tangential imposition
     * @param rData reference to element data structure
     * @return a pair of double containing the two coefficients
     */
    std::pair<const double, const double> ComputeSlipTangentialNitscheCoefficients(
        const EmbeddedElementData& rData) const;

    /**
    * This functions adds the penalty extra term level set contribution.
    * @param rLHS reference to the LHS matrix
    * @param rRHS reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddBoundaryConditionPenaltyContribution(
        MatrixType& rLHS,
        VectorType& rRHS,
        const EmbeddedElementData& rData) const;

    /**
     * This function computes the penalty coefficient for the level set BC imposition
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rData reference to element data structure
     */
    double ComputePenaltyCoefficient(
        const EmbeddedElementData& rData) const;

    /**
    * This drops the outer nodes velocity constributions in both LHS and RHS matrices.
    * @param rLHS reference to the LHS matrix
    * @param rRHS reference to the RHS vector
    * @param rData reference to element data structure
    */
    void DropOuterNodesVelocityContribution(
        MatrixType& rLHS,
        VectorType& rRHS,
        const EmbeddedElementData& rData) const;

    /**
    * This functions adds the level set strong boundary condition imposition contribution.
    * @param rLHS reference to the LHS matrix
    * @param rRHS reference to the RHS vector
    * @param rData reference to element data structure
    */
    void AddBoundaryConditionModifiedNitscheContribution(
        MatrixType& rLHS,
        VectorType& rRHS,
        const EmbeddedElementData& rData) const;

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
        EmbeddedElementData& rData,
        array_1d<double,3>& rDragForce) const;

    /**
     * @brief Calculates the location of the drag force
     * For an intersected element, this method calculates the drag force location.
     * Note that the drag force includes both the shear and the pressure contributions.
     * @param rData reference to the embedded elemental data
     * @param rDragForce reference to the computed drag force
     */
    void CalculateDragForceCenter(
        EmbeddedElementData& rData,
        array_1d<double,3>& rDragForceLocation) const;

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
    EmbeddedFluidElement& operator=(EmbeddedFluidElement const& rOther);

    /// Copy constructor.
    EmbeddedFluidElement(EmbeddedFluidElement const& rOther);

    ///@}


}; // Class EmbeddedFluidElement

namespace Internals {

template <size_t TDim, size_t TNumNodes>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator(
    const Element& rElement, const Vector& rDistance);

}

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TElementData >
inline std::istream& operator >>(std::istream& rIStream,
                                 EmbeddedFluidElement<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const EmbeddedFluidElement<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_EMBEDDED_FLUID_ELEMENT_H
