//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#ifndef KRATOS_FLUID_ELEMENT_H
#define KRATOS_FLUID_ELEMENT_H

#include "includes/checks.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/constitutive_law.h"
#include "geometries/geometry.h"

#include "includes/cfd_variables.h"
#include "custom_utilities/fluid_element_data.h"
#include "fluid_dynamics_application_variables.h"


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

// Forward decalration of auxiliary class
namespace Internals {
template <class TElementData, bool TDataKnowsAboutTimeIntegration>
class FluidElementTimeIntegrationDetail;
}

template <class TElementData>
class FluidElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FluidElement
    KRATOS_CLASS_POINTER_DEFINITION(FluidElement);

    /// Node type (default is: Node<3>)
    typedef Node<3> NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

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
    typedef MatrixRow< Matrix > ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef Kratos::Matrix ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    using ElementData = TElementData;

    static constexpr unsigned int Dim = TElementData::Dim;

    static constexpr unsigned int NumNodes = TElementData::NumNodes;

    static constexpr unsigned int BlockSize = Dim + 1;

    static constexpr unsigned int LocalSize = NumNodes * BlockSize;

    static constexpr unsigned int StrainSize = TElementData::StrainSize;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    FluidElement(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    FluidElement(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    FluidElement(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    FluidElement(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    virtual ~FluidElement();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Create a new element of this type
    /**
     * Returns a pointer to a new FluidElement element, created using given input
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
                            GeometryType::Pointer pGeom,
                            Properties::Pointer pProperties) const override;

    /// Set up the element for solution.
    /** For FluidElement, this initializes the constitutive law using the data in the element's properties.
     */
    void Initialize() override;

    /**
     * @brief CalculateLocalSystem Return empty matrices and vectors of appropriate size.
     * This element does not have a local contribution in terms of displacements, but the scheme may
     * require a proper-sized matrix, even if it is empty.
     * @param rLeftHandSideMatrix Local finite element system matrix (output)
     * @param rRightHandSideVector Local finite element residual vector (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief CalculateLeftHandSide Return an empty matrix of appropriate size.
     * This element does not have a local contribution in terms of displacements, but the scheme may
     * require a proper-sized matrix, even if it is empty.
     * @param rLeftHandSideMatrix Local finite element system matrix (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief CalculateRightHandSide Return an empty matrix of appropriate size.
     * This element does not have a local contribution in terms of displacements, but the scheme may
     * require a proper-sized matrix, even if it is empty.
     * @param rRightHandSideVector Local finite element residual vector (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief CalculateLocalVelocityContribution Calculate the local contribution in terms of velocity and pressure.
     * @param rDampMatrix Local finite element system matrix (output)
     * @param rRightHandSideVector Local finite element residual vector (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    void CalculateLocalVelocityContribution(
        MatrixType &rDampMatrix,
        VectorType &rRightHandSideVector,
        ProcessInfo &rCurrentProcessInfo) override;

    /**
     * @brief MassMatrix Calculate the local mass matrix.
     * @param rMassMatrix Local mass matrix (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    void CalculateMassMatrix(
        MatrixType &rMassMatrix,
        ProcessInfo &rCurrentProcessInfo) override;

    /**
     * @brief EquationIdVector Returns the global system rows corresponding to each local row.
     * @param rResult rResult[i] is the global index of local row i (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief GetDofList Returns a list of the element's Dofs.
     * @param rElementalDofList List of DOFs. (output)
     * @param rCurrentProcessInfo Current ProcessInfo instance. (input)
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo) override;


    /**
     * @brief GetFirstDerivativesVector Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) PRESSURE for each node.
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetFirstDerivativesVector(Vector& Values, int Step = 0) override;



    /**
     * @brief Returns ACCELERATION_X, ACCELERATION_Y, (ACCELERATION_Z,) 0 for each node
     * @param Values Vector of nodal second derivatives
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetSecondDerivativesVector(Vector& Values, int Step = 0) override;


    /**
     * @brief GetIntegrationMethod Return the integration order to be used.
     * @return Gauss Order
     */
    GeometryData::IntegrationMethod GetIntegrationMethod() const override;


    ///@}
    ///@name Access
    ///@{

    void GetValueOnIntegrationPoints(Variable<array_1d<double, 3>> const& rVariable,
                                     std::vector<array_1d<double, 3>>& rValues,
                                     ProcessInfo const& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(Variable<double> const& rVariable,
                                     std::vector<double>& rValues,
                                     ProcessInfo const& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(Variable<array_1d<double, 6>> const& rVariable,
                                     std::vector<array_1d<double, 6>>& rValues,
                                     ProcessInfo const& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(Variable<Vector> const& rVariable,
                                     std::vector<Vector>& rValues,
                                     ProcessInfo const& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(Variable<Matrix> const& rVariable,
                                     std::vector<Matrix>& rValues,
                                     ProcessInfo const& rCurrentProcessInfo) override;

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


    ///@}
    ///@name Protected Operations
    ///@{

    /// Get information from TElementData at a given point.
    /** This function serves as a wrapper so that the element does not need to
     *  know if the data is an elemental value or interpolated at the point from nodal data.
     *  @param[in] rValues The field to be read from TElementData.
     *  @param[in] rN Values of the shape functions at the desired point.
     *  @return The value evaluated at that coordinate.
     */
    virtual double GetAtCoordinate(
        const typename TElementData::NodalScalarData& rValues,
        const typename TElementData::ShapeFunctionsType& rN) const;

    /// Get information from TElementData at a given point.
    /** This function serves as a wrapper so that the element does not need to
     *  know if the data is an elemental value or interpolated at the point from nodal data.
     *  @param[in] rValues The field to be read from TElementData.
     *  @param[in] rN Values of the shape functions at the desired point.
     *  @return The value evaluated at that coordinate.
     */
    virtual array_1d<double, 3> GetAtCoordinate(
        const typename TElementData::NodalVectorData& rValues,
        const typename TElementData::ShapeFunctionsType& rN) const;

    /// Get information from TElementData at a given point.
    /** This function serves as a wrapper so that the element does not need to
     *  know if the data is an elemental value or interpolated at the point from nodal data.
     *  @param[in] rValues The field to be read from TElementData.
     *  @param[in] rN Values of the shape functions at the desired point.
     *  @return The value evaluated at that coordinate.
     */
    virtual double GetAtCoordinate(
        const double Value,
        const typename TElementData::ShapeFunctionsType& rN) const;

    virtual void CalculateMaterialResponse(TElementData& rData) const;

    /// Determine integration point weights and shape funcition derivatives from the element's geometry.
    virtual void CalculateGeometryData(Vector& rGaussWeights,
                                       Matrix& rNContainer,
                                       ShapeFunctionDerivativesArrayType& rDN_DX) const;

    /**
     * @brief Write the convective operator evaluated at this point (for each nodal funciton) to an array
     * Evaluate the convective operator for each node's shape function at an arbitrary point
     * @param rResult Output vector
     * @param rConvVel Convective velocity evaluated at the integration point
     * @param DN_DX Derivatives of shape functions evaluated at the integration point
     */
    void ConvectionOperator(Vector& rResult,
                            const array_1d<double,3>& rConvVel,
                            const ShapeFunctionDerivativesType& DN_DX) const;

    virtual void AddTimeIntegratedSystem(
        TElementData& rData,
        MatrixType& rLHS,
        VectorType& rRHS);

    virtual void AddTimeIntegratedLHS(
        TElementData& rData,
        MatrixType& rLHS);

    virtual void AddTimeIntegratedRHS(
        TElementData& rData,
        VectorType& rRHS);

    virtual void AddVelocitySystem(
        TElementData& rData,
        MatrixType& rLocalLHS,
        VectorType& rLocalRHS);

    virtual void AddMassLHS(
        TElementData& rData,
        MatrixType& rMassMatrix);

    /**
     * @brief Adds the boundary traction component along a cut plane for embedded formulations.
     * This method adds the boundary traction component to the LHS and RHS arrays.
     * Such boundary integral must be implemented in all the fluid dynamics elements
     * deriving from this one in accordance to the formulation used. This method is
     * intended to be called from the derived elements to add the contribution of the
     * tractions on the elemental cuts to enforce equilibrium. This means that what we
     * call external traction is nothing but minus the base formulation boundary term.
     * @param rData Element data structure
     * @param rUnitNormal Outwards unit normal vector for the cut plane
     * @param rLHS Reference to the Left Hand Side matrix
     * @param rRHS Reference to the Right Hand Side vector
     */
    virtual void AddBoundaryTraction(
        TElementData& rData,
        const Vector& rUnitNormal,
        MatrixType& rLHS,
        VectorType& rRHS);

    void GetCurrentValuesVector(
        const TElementData& rData,
        array_1d<double,LocalSize>& rValues) const;

    ///@}
    ///@name Protected  Access
    ///@{

    const ConstitutiveLaw::Pointer GetConstitutiveLaw() const;

    ConstitutiveLaw::Pointer GetConstitutiveLaw();

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

    //// Constitutive relation for the element
    ConstitutiveLaw::Pointer mpConstitutiveLaw = nullptr;

    ///@}
    ///@name Friends
    ///@{

    friend class Internals::FluidElementTimeIntegrationDetail<TElementData, TElementData::ElementManagesTimeIntegration>;

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
    FluidElement& operator=(FluidElement const& rOther);

    /// Copy constructor.
    FluidElement(FluidElement const& rOther);

    ///@}


}; // Class FluidElement

namespace Internals {

template< class TElementData, std::size_t TDim >
struct StrainRateSpecialization {
/// Compute the strain rate vector in Voigt notation, to use as input for the fluid constitutive law.
/*  @param[out] rStrainRate The strain rate tensor (symmetric gradient of velocity) in Voigt notation.
 *  @param[in] rVelocities Matrix of nodal velocities, as provided by TElementData.
 *  @param[in] rDNDX Matrix of shape function gradients on the integration point, as provided by TElementData.
 *  @see ConstitutiveLaw.
 */
static void Calculate(
    Vector& rStrainRate,
    const typename TElementData::NodalVectorData& rVelocities,
    const typename TElementData::ShapeDerivativesType& rDNDX);
};

template< class TElementData >
struct StrainRateSpecialization<TElementData,2> {
static void Calculate(
    Vector& rStrainRate,
    const typename TElementData::NodalVectorData& rVelocities,
    const typename TElementData::ShapeDerivativesType& rDNDX);
};

template< class TElementData >
struct StrainRateSpecialization<TElementData,3> {
static void Calculate(
    Vector& rStrainRate,
    const typename TElementData::NodalVectorData& rVelocities,
    const typename TElementData::ShapeDerivativesType& rDNDX);
};

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
                                 FluidElement<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const FluidElement<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_FLUID_ELEMENT_H
